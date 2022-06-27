#!/usr/bin/env python

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
import os


class TestNWEuropeanShelf(unittest.TestCase):
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
        #Filename of vector (shapefile) containing all the boundaries
        # for this test-suite. Annotated with PhysID.
        self.boundaries_filename = os.path.join(
            self.this_path, 'test_NWEuropeanShelf_lines.shp')
        #Filename of vector (shapefile) containing all the polygons
        # for this test-suite.
        self.polygons_filename = os.path.join(
            self.this_path, 'test_NWEuropeanShelf_polygons.shp')
        shapefile_set=[".shp", ".shx", ".prj", ".qpj", ".dbf", ".cpg"]
        for extension in shapefile_set:
            self.output_files.append(
                    os.path.splitext(self.polygons_filename)[0] +
                    extension)
        #Filename of gmsh geometry definition.
        self.mesh_generator_geometry = os.path.join(
                self.this_path, 'test_NWEuropeanShelf.geo')
        self.output_files.append(self.mesh_generator_geometry)
        #Mesh filename
        self.mesh_filename = os.path.join(
                self.this_path, 'test_NWEuropeanShelf.msh')
        self.output_files.append(self.mesh_filename)


    def remove_file(self, filename):
        '''Remove files, whithout failing for files that do not exist.

        This method uses os.remove to remove files, but captures the exception
        thrown in case the files do not exist. Exceptions due to other reasons
        are still raised.
        '''
        import qmesh
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
        except:
            raise
        else:
            qmesh.LOG.info('Removed file '+filename)


    def tearDown(self):
        '''Method taking down the test fixture, tiding up.

        Test-cases typically generate files that can obscure working. This method
        will remove all files identified in the setUp method as output files. Note that
        the output files are ignored by the git set-up.
        '''
        for filename in self.output_files:
            self.remove_file(filename)


    def test_NWEuropeanShelf_UTM30(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh
        except ImportError:
            os.sys.path.append(self.project_root_path)
            import qmesh
        #Try initialising qgis API
        try:
            qmesh.initialise()
            self.assertTrue(True)
        except AssertionError:
            self.assertTrue(False)
 
        #Try reading-in prepared boundaries (lines)
        try:
            boundaries = qmesh.vector.Shapes()
            boundaries.fromFile(self.boundaries_filename)
        except:
            self.assertTrue(False)
 
        #Try assembling the region (polygon) from lines
        try:
            loops = qmesh.vector.identifyLoops(boundaries,fixOpenLoops=True)
            poygons = qmesh.vector.identifyPolygons(loops)
            poygons.writeInfo()
            poygons.writeFile(self.polygons_filename)
            self.assertTrue(True)
        except AssertionError:
            self.assertTrue(False)

        #Convert to Gmsh geometry
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, poygons)
            domain.setTargetCoordRefSystem('EPSG:32630')
        except AssertionError:
            self.assertTrue(False)
        #Try meshing. Note, qmesh.gmshDomain() defaults to del2d meshing
        # algo, but in this case it does not give good results.
        try:
            domain.gmsh(geoFilename=self.mesh_generator_geometry,
                        mshFilename=self.mesh_filename,
                        )
        except AssertionError:
            self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
