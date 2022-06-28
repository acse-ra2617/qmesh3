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

import unittest
import sys
import os
import errno
import glob

class Test_Shapes(unittest.TestCase):
    '''Todo: add test documentation as class docstring '''
    def setUp(self):
        '''Method setting-up test fixture.

        This method is not aimed to be invoked by the user. Various resources needed and
        produced by the test are defined here:
        1) The path of the testing function and the project root path (qmesh).
        2) The filenames, with absolute paths, of files needed by the test-cases.
        3) The filenames, with absolute paths, of files produced by the test-cases.
        '''
        #Define list of output files, empty for now
        self.output_files = []
        #Construct absolute path to this test file, without the test file-name
        self.this_path = os.path.dirname(os.path.realpath(__file__))
        #Construct absolute path to the test-type directory
        self.project_test_type_path = os.path.split(self.this_path)[0]
        #Construct absolute path to the tests directory
        self.project_tests_path = os.path.split(self.project_test_type_path)[0]
        #Construct absolute path to the project root
        self.project_root_path = os.path.split(self.project_tests_path)[0]
        #Define filenames used in this test suite.
        #Filename of vector (shapefile) containing multi-lines,
        # to be decomposed into lines by this test
        self.multiLines_filename = os.path.join(self.this_path,'multiLines.shp')
        #
        self.multiPolygons_filename = os.path.join(self.this_path,'multiPolygons.shp')
        #
        self.lines_filename = os.path.join(self.this_path,'lines')
        self.output_files.extend(glob.glob(self.lines_filename+'.*'))
        #
        self.polygons_filename = os.path.join(self.this_path,'polygons.shp')
        self.output_files.extend(glob.glob(self.polygons_filename+'.*'))

    def remove_file(self, filename):
        '''Remove files, whithout failing for files that do not exist.

        This method uses os.remove to remove files, but captures the exception
        thrown in case the files do not exist. Exceptions due to other reasons
        are still raised.
        '''
        try:
            os.remove(filename)
        except OSError as error:
            if error.errno != errno.ENOENT:
                raise

    def tearDown(self):
        '''Method taking down the test fixture, tiding up.

        Test-cases typically generate files that can obscure working. This method
        will remove all files identified in the setUp method as output files. Note that
        the output files are ignored by the git set-up.
        '''
        for filename in self.output_files:
            self.remove_file(filename)

    def test_decomposeMultiLines(self):
        """ Test successful conversion of multilines into lines. """
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        qmesh3.initialise()
        lines = qmesh3.vector.Shapes()
        lines.fromFile(self.multiLines_filename)
        lines.decomposeMultiFeatures()
        lines.writeFile(self.lines_filename)
        
    def test_decomposeMultiPolygons(self):
        """ Test successful conversion of multipolygons into polygons. """
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        qmesh3.initialise()
        polygons = qmesh3.vector.Shapes()
        polygons.fromFile(self.multiPolygons_filename)
        polygons.decomposeMultiFeatures()
        polygons.writeFile(self.polygons_filename)

suite = unittest.TestLoader().loadTestsFromTestCase(Test_Shapes)

if __name__ == '__main__':
    unittest.main()
