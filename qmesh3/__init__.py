#    Copyright (C) 2013 Alexandros Avdis and others. See the AUTHORS file for a full list of copyright holders.
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
from  qgis.core import QgsApplication

import qmesh3.vector
import qmesh3.raster
import qmesh3.mesh
import qmesh3.publish
import qmesh3.lib
from qmesh3.config import *

import pkg_resources
#Set the version attribute
try:
    __packaged_distro__ = pkg_resources.get_distribution('qmesh3')
    __version__ = __packaged_distro__.version
except (pkg_resources.DistributionNotFound, AttributeError, NameError):
    __version__ = "local"
    LOG.warning('Could not find qmesh version information. Provenance information is incomplete.')
#Set the git-sha-key attribute
try:
    __packaged_distro__.has_metadata('GIT_SHA_KEY')
    __git_sha_key__ = __packaged_distro__.get_metadata('GIT_SHA_KEY')
except (AttributeError, FileNotFoundError, NameError):
    __git_sha_key__ = "local"
    LOG.warning('Could not find qmesh origin git SHA key. Provenance information is incomplete.')
#Set the license attribute
try:
    __packaged_distro__.has_metadata('LICENSE')
    __license__ = __packaged_distro__.get_metadata('LICENSE')
except (AttributeError, FileNotFoundError, NameError):
    __license__ = ""
    LOG.warning('Could not find complete license statement. Provenance information is incomplete.')
    LOG.warning('Qmesh is distributed under GPLv3. Please observe the license at .')
#Set the authors attribute
try:
    __packaged_distro__.has_metadata('AUTHORS.md')
    __authors__ = __packaged_distro__.get_metadata('AUTHORS.md')
except (AttributeError, FileNotFoundError, NameError):
    __authors__ = ""
    LOG.warning('Could not find qmesh authors information. Provenance information is incomplete.')
# Set the gmsh-path attribute
try:
    __packaged_distro__.has_metadata('GMSH_BIN_PATH')
    __gmsh_bin_path__ = __packaged_distro__.get_metadata('GMSH_BIN_PATH')
except (AttributeError, FileNotFoundError, NameError):
    __gmsh_bin_path__ = "/usr/bin"
    LOG.warning('Could not find gmsh binary path specification. Mesh generation might fail.')
# Set the qgis-install-path attribute
try:
    __packaged_distro__.has_metadata('QGIS_PATH')
    __qgis_path__ = __packaged_distro__.get_metadata('QGIS_PATH')
except (AttributeError, FileNotFoundError, NameError):
    __qgis_path__ = "/usr/local/"
    LOG.warning('Could not find qgis path specification. qgis initialisation might fail.')
try:
    del __packaged_distro__
except NameError:
    pass

# Set-up environment and initialise qgis, without graphics.
# Suppress verbose debugging Qt messages
os.environ['QT_LOGGING_RULES']="qt5ct.debug=false"
# Make sure qmesh can run in environments without graphics
os.environ["QT_QPA_PLATFORM"] = "offscreen"
# Initialise qgis
qgs = QgsApplication([], True, None)
qgis_install_path='/usr'
qgs.setPrefixPath(__qgis_path__, True)
qgs.initQgis()
