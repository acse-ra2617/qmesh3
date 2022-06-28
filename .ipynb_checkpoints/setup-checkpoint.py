#!/usr/bin/env python
#
#    Copyright (C) 2013 Alexandros Avdis and others.
#    See the AUTHORS.md file for a full list of copyright holders.
#
#    This file is part of QMesh.
#
#    QMesh is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    QMesh is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with QMesh.  If not, see <http://www.gnu.org/licenses/>.

import os
import io
from setuptools import setup
from setuptools.command.install import install as setuptools_install
from setuptools.command.develop import develop as setuptools_develop

# Suppress verbose debugging Qt messages
os.environ['QT_LOGGING_RULES'] = "qt5ct.debug=false"
# Make sure qmesh can run in environments without graphics
os.environ["QT_QPA_PLATFORM"] = "offscreen"

def main():
    current_directory = os.path.dirname(__file__)
    # Read version from file. The variable storing the
    # version identifier is passed into a standardised argument of 
    # setuptools.setup
    ##version_path = os.path.join(current_directory, 'VERSION')
    ##with io.open(version_path, encoding='utf-8') as version_file:
        ##qmesh_version = version_file.readline().strip()
    # Read long description from file. The variable storing the
    # long description is passed into a standardised argument of 
    # setuptools.setup
    readme_path = os.path.join(current_directory, 'README.rst')
    with io.open(readme_path, encoding='utf-8') as readme_file:
        long_description = readme_file.read()

        
    setup(
          name='qmesh3',
          version='1.0.3',
          description = "Finite Element meshes from GIS data.",
          long_description = long_description,
          author = "The QMesh Development Team.",
          author_email = "develop@qmesh.org",
          url = "http://www.qmesh.org",
          download_url = 'https://github.com/acse-ra2617/qmesh3/releases/tag/1.0.3',
          packages = ['qmesh3',
                      'qmesh3.lib',
                      'qmesh3.vector',
                      'qmesh3.mesh',
                      'qmesh3.raster',
                      'qmesh3.publish',
                     ],
          package_dir = {
              'qmesh3': 'qmesh3',
              'qmesh3.lib':'qmesh3/lib',
              'qmesh3.vector':'qmesh3/vector',
              'qmesh3.meshg':'qmesh3/mesh',
              'qmesh3.raster':'qmesh3/raster',
              'qmesh3.publish':'qmesh3/publish',
              },
          provides=['qmesh3'],
          install_requires=['setuptools-qmesh', 'GFD_basisChangeTools'],
          setup_requires=['setuptools>=35.0.0', 'setuptools_qmesh'],
          extras_require={'RDM':['pyrdm']},
          entry_points={
            "distutils.commands":[
                "check_qgis = setuptools_qmesh.command.check_qgis:CheckQgis",
                "check_gmsh = setuptools_qmesh.command.check_gmsh:CheckGmsh"],
            "distutils.setup_keywords": [
                "qgis_path            = setuptools_qmesh.dist:assert_path",
                "gmsh_bin_path        = setuptools_qmesh.dist:assert_path",
                "include_git_sha_key  = setuptools.dist:assert_bool",
                "include_full_license = setuptools.dist:assert_bool",
                "include_author_ids   = setuptools.dist:assert_bool"],
            "egg_info.writers": [
                "QGIS_PATH     = setuptools_qmesh.command.egg_info:write_qgis_path",
                "GMSH_BIN_PATH = setuptools_qmesh.command.egg_info:write_gmsh_bin_path",
                "GIT_SHA_KEY   = setuptools_qmesh.command.egg_info:write_git_sha_key",
                "LICENSE       = setuptools_qmesh.command.egg_info:write_full_license",
                "AUTHORS.md    = setuptools_qmesh.command.egg_info:write_author_ids"],
            },
          include_git_sha_key=True,
          include_full_license=True,
          include_author_ids=True,
          license='GPLv3',
          test_suite = "tests",
          keywords = ['GIS', 'mesh generation'],
          cmdclass={'install': Install,
                    'develop': Develop}
        )

class Install(setuptools_install):
    '''Class invoking installation procedure.

    This class is overloading the setuptools install class. The
    run method of the parent class is run.
    '''
    def run(self):
        from setuptools_qmesh.command.check_gmsh import CheckGmsh
        from setuptools_qmesh.command.check_qgis import CheckQgis
        #Check gmsh is installed
        gmsh_checker = CheckGmsh(self.distribution)
        gmsh_checker.initialize_options()
        gmsh_checker.finalize_options()
        gmsh_checker.run()
        #Check qgis is installed
        qgis_checker = CheckQgis(self.distribution)
        qgis_checker.initialize_options()
        qgis_checker.finalize_options()
        qgis_checker.run()
        #Install qmesh
        setuptools_install.run(self)

class Develop(setuptools_develop):
    '''Class invoking installation-in-development-mode procedure.

    Installation in development mode does not copy files to
    standardised locations. Files in standardised locations are
    still created, but they point to files in the directory where
    the source code resides. Thus, the source code in the working
    directory is run as the installed code. Modifications are
    immediately propagated to the installation, and development
    is facilitated. This class is overloading the setuptools
    develop class. The run method of the parent class is run.
    '''
    def run(self):
        from setuptools_qmesh.command.check_gmsh import CheckGmsh
        from setuptools_qmesh.command.check_qgis import CheckQgis
        #Check gmsh is installed
        gmsh_checker = CheckGmsh(self.distribution)
        gmsh_checker.initialize_options()
        gmsh_checker.finalize_options()
        gmsh_checker.run()
        #Check qgis is installed
        qgis_checker = CheckQgis(self.distribution)
        qgis_checker.initialize_options()
        qgis_checker.finalize_options()
        qgis_checker.run()
        #Install qmesh
        setuptools_develop.run(self)

if __name__=='__main__':
    main()
