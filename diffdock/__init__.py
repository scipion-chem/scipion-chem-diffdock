# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using ConPLex models for virtual screening
"""

# General imports
import os, subprocess, json

# Scipion em imports
import pwem
from scipion.install.funcs import InstallHelper

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'mit_logo.png'

class Plugin(pwchemPlugin):
	"""
	"""
	_dfdHome = os.path.join(pwem.Config.EM_ROOT, DIFFDOCK_DIC['name'] + '-' + DIFFDOCK_DIC['version'])

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(DIFFDOCK_DIC['home'], cls._dfdHome)

	@classmethod
	def defineBinaries(cls, env):
		"""
        This function defines the binaries for each package.
        """
		cls.addDiffDockPackage(env)

	@classmethod
	def addDiffDockPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		# Instantiating the install helper
		installer = InstallHelper(DIFFDOCK_DIC['name'], packageHome=cls.getVar(DIFFDOCK_DIC['home']),
															packageVersion=DIFFDOCK_DIC['version'])

		# Installing package
		installer.getCloneCommand(cls.getDiffDockGithub(), targeName='DIFFDOCK_CLONED') \
			.getCondaEnvCommand(pythonVersion='3.9', requirementsFile=False) \
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && '
		                f'pip install torch==1.13.1+cu117 '
		                f'--extra-index-url https://download.pytorch.org/whl/cu117', 'PYTORCH_INSTALLED') \
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && '
						f'conda install -y -c conda-forge prody==2.2.0 && '
		                f'pip install torch-cluster==1.6.0+pt113cu117 torch-sparse==0.6.16+pt113cu117 '
		                f'torch-scatter==2.1.0+pt113cu117 torch-spline-conv==1.2.1+pt113cu117 '
		                f'torch-geometric==2.2.0 '
		                f'--find-links https://pytorch-geometric.com/whl/torch-1.13.1+cu117.html', 'DIFFDOCK_INSTALLED') \
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && '
		                f'pip install e3nn==0.5.1 fair-esm==2.0.0 networkx==2.8.4 pandas==1.5.1 '
		                f'pybind11==2.11.1 pytorch-lightning==1.9.5 rdkit==2022.03.3 '
		                f'scikit-learn==1.1.0 torchmetrics==0.11.0 dllogger@git+https://github.com/NVIDIA/dllogger.git  '
		                f'biopython PyYAML scipy spyrmsd biopandas', 'ESM_INSTALLED') \
			.addPackage(env, ['git', 'conda', 'pip'], default=default)


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def getPackageDir(cls, path=''):
		return os.path.abspath(os.path.join(cls.getVar(DIFFDOCK_DIC['home']), path))

	@classmethod
	def getDiffDockGithub(cls):
		return 'https://github.com/gcorso/DiffDock.git'

	@classmethod
	def getPluginHome(cls, path=""):
		import diffdock
		fnDir = os.path.split(diffdock.__file__)[0]
		return os.path.join(fnDir, path)