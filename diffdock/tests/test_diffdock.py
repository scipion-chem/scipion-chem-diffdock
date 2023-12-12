# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb, ProtSetFilter

from pwchem.protocols import ProtChemImportSmallMolecules

from ..protocols import ProtDiffDockDocking

class TestDiffDock(BaseTest):
  @classmethod
  def setUpClass(cls):
    cls.ds = DataSet.getDataSet('model_building_tutorial')
    cls.dsLig = DataSet.getDataSet("smallMolecules")
    setupTestProject(cls)

    cls._runImportPDB()
    cls._runImportSmallMols()
    cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
    cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

  @classmethod
  def _runImportPDB(cls):
    cls.protImportPDB = cls.newProtocol(
      ProtImportPdb,
      inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
    cls.proj.launchProtocol(cls.protImportPDB, wait=False)

  @classmethod
  def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
          ProtChemImportSmallMolecules,
          filesPath=cls.dsLig.getFile('mol2'))
      cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

  def _runSetFilter(self, inProt, number, property):
    self.protFilter = self.newProtocol(
      ProtSetFilter,
      operation=ProtSetFilter.CHOICE_RANKED,
      threshold=number, rankingField=property)
    self.protFilter.inputSet.set(inProt)
    self.protFilter.inputSet.setExtended('outputSmallMolecules')

    self.proj.launchProtocol(self.protFilter, wait=False)
    return self.protFilter

  def _runDiffDock(self, recProt, ligProt):
    protDiffDock = self.newProtocol(
    ProtDiffDockDocking,
    inputAtomStruct=recProt.outputPdb,
    inputSmallMols=ligProt.outputSmallMolecules,
    nSamples=5, inferSteps=10)
    self.proj.launchProtocol(protDiffDock, wait=False)

    return protDiffDock

  def test(self):
    protSmallFilter = self._runSetFilter(inProt=self.protImportSmallMols, number=2, property='smallMoleculeFile')
    self._waitOutput(protSmallFilter, 'outputSmallMolecules', sleepTime=10)

    print('Docking with DiffDock in the whole protein')
    protDiffDock = self._runDiffDock(self.protImportPDB, protSmallFilter)

    self._waitOutput(protDiffDock, 'outputSmallMolecules', sleepTime=10)
    self.assertIsNotNone(getattr(protDiffDock, 'outputSmallMolecules', None))

