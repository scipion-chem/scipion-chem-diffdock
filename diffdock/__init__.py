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

		ymlFile = "environment.yml"
		splitScript = cls.getScriptsDir('splitPipBlocks.py')

		# Installing package
		installer.getCloneCommand(cls.getDiffDockGithub(), binaryFolderName='.', targeName='DIFFDOCK_CLONED') \
			.addCommand(f"sed -i 's/name: diffdock/name: {cls.getEnvName(DIFFDOCK_DIC)}/g' {ymlFile} && "
									f'conda env create -f "$(python {splitScript} {ymlFile} True)" && '
									f'python {splitScript} {ymlFile} False | while read -r f; do conda env update -f "$f"; done',
									'DIFFDOCK_INSTALLED')\
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
