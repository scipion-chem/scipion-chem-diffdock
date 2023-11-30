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
		installer.getCloneCommand('https://github.com/gcorso/DiffDock.git', targeName='DIFFDOCK_CLONED') \
			.getCondaEnvCommand(pythonVersion='3.9', requirementsFile=False) \
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && conda install -y pytorch==1.11.0 pytorch-cuda=11.7 -c pytorch -c nvidia', 'PYTORCH_INSTALLED')\
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && '
									f'pip install --jobs=6 torch-scatter torch-sparse torch-cluster torch-spline-conv torch-geometric==2.0.4 '
									f'-f https://data.pyg.org/whl/torch-1.11.0+cu117.html && '
									f'python -m pip install --jobs=6 PyYAML scipy "networkx[default]" '
									f'biopython rdkit-pypi e3nn spyrmsd pandas biopandas', 'DIFFDOCK_INSTALLED') \
			.addCommand(f'{cls.getEnvActivationCommand(DIFFDOCK_DIC)} && pip install "fair-esm[esmfold]" && '
									f'pip install "dllogger @ git+https://github.com/NVIDIA/dllogger.git" && '
									f'pip install "openfold @ git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307"', 'ESM_INSTALLED') \
			.addPackage(env, ['git', 'conda', 'pip'], default=default)


	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runScript(cls, protocol, scriptName, args, envDict, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDict), 'python', scriptName)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)
