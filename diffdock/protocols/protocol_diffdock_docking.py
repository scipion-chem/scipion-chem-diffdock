# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pwem.convert.atom_struct import toPdb

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseName, pdbqt2other

from .. import Plugin as diffdockPlugin
from ..constants import DIFFDOCK_DIC

class ProtDiffDockDocking(EMProtocol):
  """Run a prediction using a DiffDock trained model over a proteins and a set of ligands"""
  _label = 'diffdock docking'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                    label='Input atomic structure: ',
                    help="The atomic structure to use as receptor in the docking")
    iGroup.addParam('inputSmallMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                    label='Input small molecules: ',
                    help='Set of small molecules to input the model for predicting their interactions')

    mGroup = form.addGroup('Model', expertLevel=params.LEVEL_ADVANCED)
    mGroup.addParam('scoreModel', params.PathParam, label='Score model (pt): ', default='',
                    help="Select the model file you want to use for scoring. "
                         "\nIf None, the default behaviour will be used.")
    mGroup.addParam('confidenceModel', params.PathParam, label='Confidence model (pt): ', default='',
                    help="Select the model file you want to use for confidence calculation. "
                         "\nIf None, the default behaviour will be used.")

    pGroup = form.addGroup('Prediction')
    pGroup.addParam('nSamples', params.IntParam, label='Number of positions: ', default=20,
                    help='Number of positions to generate')
    pGroup.addParam('inferSteps', params.IntParam, label='Inference steps: ', default=20,
                   help='Number of denoising steps')

    pGroup.addParam('batchSize', params.IntParam, label='Batch size: ', default=10, expertLevel=params.LEVEL_ADVANCED,
                    help='Batch size for the model')
    pGroup.addParam('finalDenoise', params.BooleanParam, label='Final step denoise: ', default=False,
                    expertLevel=params.LEVEL_ADVANCED,
                    help='Whether to use no noise in the final step of the reverse diffusion')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.convertStep)
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.createOutputStep)


  def convertStep(self):
    smiDir = self.getInputSMIDir()
    if not os.path.exists(smiDir):
      os.makedirs(smiDir)

    molDir = self.copyInputMolsInDir()
    args = f' --multiFiles -iD "{molDir}" --pattern "*" -of smi --outputDir "{smiDir}"'
    pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=smiDir)

    inASFile = self.inputAtomStruct.get().getFileName()
    outASFile = os.path.abspath(self._getTmpPath(getBaseName(inASFile) + '.pdb'))
    if inASFile.endswith('.pdbqt'):
      pdbqt2other(self, inASFile, outASFile)
    else:
      toPdb(inASFile, outASFile)

  def predictStep(self):
    csvFile = self.buildCSVFile()
    outDir = os.path.abspath(self._getExtraPath())

    program = f'{pwchemPlugin.getEnvActivationCommand(DIFFDOCK_DIC)} && python -m inference '
    args = f'--protein_ligand_csv {csvFile} --out_dir {outDir} '
    args += f'--inference_steps {self.inferSteps.get()} --samples_per_complex {self.nSamples.get()} ' \
            f'--batch_size {self.batchSize.get()} '
    if not self.finalDenoise.get():
      args += '--no_final_step_noise '

    scoreModelDir = os.path.dirname(self.scoreModel.get()) if self.scoreModel.get() else None
    confModelDir = os.path.dirname(self.confidenceModel.get()) if self.confidenceModel.get() else None
    if scoreModelDir:
      args += f'--model_dir {scoreModelDir} '
    if confModelDir:
      args += f'--confidence_model_dir {confModelDir} '

    self.runJob(program, args, cwd=diffdockPlugin.getPackageDir())

  def createOutputStep(self):
    outDir = self._getPath('outputLigands')
    if not os.path.exists(outDir):
      os.mkdir(outDir)
    outputSet = SetOfSmallMolecules().create(outputPath=outDir)

    outDic = self.parseOutputDocks()
    for smallMol in self.inputSmallMols.get():
      molName = getBaseName(smallMol.getFileName())
      for outFile in outDic[molName]:
        conf = outFile.split('_confidence')[-1].split('.sdf')[0]
        posId = outFile.split('/rank')[-1].split('_')[0]

        newSmallMol = SmallMolecule()
        newSmallMol.copy(smallMol, copyId=False)
        newSmallMol._energy = pwobj.Float(conf)
        newSmallMol.poseFile.set(outFile)
        newSmallMol.setPoseId(posId)
        newSmallMol.gridId.set(1)
        newSmallMol.setMolClass('DiffDock')
        newSmallMol.setDockId(self.getObjId())

        outputSet.append(newSmallMol)


    outputSet.proteinFile.set(self.inputAtomStruct.get().getFileName())
    outputSet.setDocked(True)
    self._defineOutputs(outputSmallMolecules=outputSet)


  ###########################################################

  def parseOutputDocks(self):
    outDirs = []
    for oDir in os.listdir(self._getExtraPath()):
      if os.path.isdir(self._getExtraPath(oDir)) and oDir != 'inputSMI':
        outDirs.append(oDir)

    outDic = {}
    for oDir in outDirs:
      outDic[oDir] = []
      for outFile in os.listdir(self._getExtraPath(oDir)):
        if '_confidence' in outFile:
          outDic[oDir].append(os.path.join(self._getExtraPath(oDir), outFile))

    return outDic

  def copyInputMolsInDir(self):
    oDir = os.path.abspath(self._getTmpPath('inMols'))
    if not os.path.exists(oDir):
      os.makedirs(oDir)

    for mol in self.inputSmallMols.get():
      os.link(mol.getFileName(), os.path.join(oDir, os.path.split(mol.getFileName())[-1]))
    return oDir

  def getInputSMIDir(self):
    return os.path.abspath(self._getExtraPath('inputSMI'))

  def getInputSMIs(self):
    smisDic = {}
    iDir = self.getInputSMIDir()
    for file in os.listdir(iDir):
      with open(os.path.join(iDir, file)) as f:
        title, smi = f.readline().split()
        smisDic[title] = smi.strip()
    return smisDic
  
  def getInputASFile(self):
    iASFile = None
    for file in os.listdir(self._getTmpPath()):
      if file.endswith('.pdb'):
        iASFile = os.path.abspath(self._getTmpPath(file))
    return iASFile

  def getInputCSV(self):
    return os.path.abspath(self._getExtraPath('inputPairs.csv'))

  def buildCSVFile(self):
    csvFile = self.getInputCSV()
    smiDic = self.getInputSMIs()
    iASFile = self.getInputASFile()

    with open(csvFile, 'w') as f:
      f.write('complex_name,protein_path,ligand_description,protein_sequence\n')
      for smi, title in smiDic.items():
        f.write(f'{title},{iASFile},{smi},\n')
    return csvFile
