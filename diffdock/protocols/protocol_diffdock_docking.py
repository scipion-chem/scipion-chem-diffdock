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

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.utils import getBaseName, pdbqt2other

from .. import Plugin as diffdockPlugin
from ..constants import DIFFDOCK_DIC

class ProtDiffDockDocking(EMProtocol):
  """Run a prediction using a ConPLex trained model over a set of proteins and ligands"""
  _label = 'diffdock virtual screening'

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

    pGroup = form.addGroup('Prediction')
    pGroup.addParam('nSamples', params.IntParam, label='Number of positions: ', default=20,
                    help='Number of positions to generate')
    pGroup.addParam('inferSteps', params.IntParam, label='Inference steps: ', default=20,
                   help='Number of denoising steps')

    pGroup.addParam('batchSize', params.IntParam, label='Batch size: ', default=10, expertLevel=params.LEVEL_ADVANCED,
                    help='Batch size for the model')
    pGroup.addParam('actualSteps', params.IntParam, label='Actual steps: ', default=18,
                    expertLevel=params.LEVEL_ADVANCED, help='Number of denoising steps that are actually performed')
    pGroup.addParam('finalDenosise', params.BooleanParam, label='Final step denoise: ', default=False,
                    expertLevel=params.LEVEL_ADVANCED,
                    help='Whether to use no noise in the final step of the reverse diffusion')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.convertStep)
    self._insertFunctionStep(self.predictStep)
    # todo: how to express the output
    # self._insertFunctionStep('createOutputStep')


  def convertStep(self):
    outDir = self.getInputMolsDir()
    for mol in self.inputSmallMols.get():
      fnSmall = os.path.abspath(mol.getFileName())
      if fnSmall.endswith('.pdbqt'):
        args = f' -i "{fnSmall}" -of sdf --outputDir "{outDir}" --outputName {getBaseName(fnSmall)}'
        pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir, popen=True)
      else:
        os.link(fnSmall, os.path.join(outDir, os.path.split(fnSmall)[-1]))

    inASFile = self.inputAtomStruct.get().getFileName()
    outASFile = os.path.abspath(self._getTmpPath(getBaseName(inASFile) + '.pdb'))
    if inASFile.endswith('.pdbqt'):
      pdbqt2other(self, inASFile, outASFile)

  def predictStep(self):
    csvFile = self.buildCSVFile()
    outDir = os.path.abspath(self._getExtraPath())
    scoreModelPath = os.path.join(diffdockPlugin.getModelsDir('paper_score_model/'))
    confModelPath = os.path.join(diffdockPlugin.getModelsDir('paper_score_model/'))

    program = f'{pwchemPlugin.getEnvActivationCommand(DIFFDOCK_DIC)} && python -m inference '
    args = f'--protein_ligand_csv {csvFile} --out_dir {outDir} '
    args += f'--model_dir {scoreModelPath} --confidence_model_dir {confModelPath} '
    args += f'--inference_steps {self.inferSteps.get()} --samples_per_complex {self.nSamples.get()} ' \
            f'--batch_size {self.batchSize.get()} --actual_steps {self.actualSteps.get()} '
    if not self.finalDenoise.get():
      args += '--no_final_step_noise '

    self.runJob(program, args, cwd=self._getExtraPath())


  def getInputMolsDir(self):
    return os.path.abspath(self._getTmpPath('inMols'))
  
  def getInputMolFiles(self):
    imolFiles = []
    iDir = self.getInputMolsDir()
    for file in os.listdir(iDir):
      imolFiles.append(os.path.join(iDir, file))
    return imolFiles

  def getInputASFile(self):
    iASFile = None
    for file in os.listdir(self._getTmpPath()):
      if file.endswith('.pdb'):
        iASFile = self._getTmpPath(file)
    return iASFile

  def getInputCSV(self):
    return os.path.abspath(self._getExtraPath('inputPairs.csv'))

  def buildCSVFile(self):
    csvFile = self.getInputCSV()
    iMolFiles = self.getInputMolFiles()
    iASFile = self.getInputASFile()

    with open(csvFile, 'w') as f:
      f.write('complex_name,protein_path,ligand_description,protein_sequence\n')
      for molFile in iMolFiles:
        cName = getBaseName(molFile)
        f.write(f'{cName},{iASFile},{molFile},\n')
    return csvFile
