#    Copyright (C) 2013 Alexandros Avdis and others.
#    See the AUTHORS file for a full list of copyright holders.
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

class TestGGeometryPoints(unittest.TestCase):

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

    def test_add_point(self):
        '''Testing the point-addition method in qmesh.mesh.Geometry class
        '''
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        # This should fail, the point is only 2D
        try:
            myGeometry.addPoint([0,1])
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    # again, should fail as we're passing dodgy data
    def test_add_incorrect_point(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        try:
            myGeometry.addPoint(1)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    # again, should fail as we're passing dodgy data
    def test_add_point_as_string(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        try:
            myGeometry.addPoint("hhh")
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_add_existing_point(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        myGeometry.addPoint([0,1,1])
        self.assertTrue(myGeometry.points == {0:[0,1,1], 1:[0,1,1]})

    def test_remove_existing_point(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        myGeometry.removePoint(0)
        self.assertTrue(myGeometry.points == {})

    def test_remove_nonexisting_point(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        try:
            myGeometry.removePoint(1)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)
        
    def test_get_point(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        myGeometry.addPoint([0,1,2])
        self.assertTrue(myGeometry.getPoint(0) == [0,1,1])
        self.assertTrue(myGeometry.getPoint(1) == [0,1,2])
        try:
            self.assertTrue(myGeometry.getPoint(2))
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_get_point_odd_input(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        myGeometry.addPoint([0,1,2])
        try:
            myGeometry.getPoint("ibiv")
        except TypeError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_add_points_from_dict(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        points = {0: [0,1,0], 1:[1,1,2], 2:[3,54,1]}
        myGeometry.addPointsFromDictionary(points)
        self.assertTrue(myGeometry.getPoint(0) == [0,1,0])

    def test_add_points_from_dict_existing(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        points = {0: [0,1,0], 1:[1,1,2], 2:[3,54,1]}
        myGeometry.addPointsFromDictionary(points)
        self.assertTrue(myGeometry.getPoint(0) == [0,1,1])
        self.assertTrue(myGeometry.pointsCount() == 4)

    def test_add_points_existing(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        myGeometry = qmesh3.mesh.Geometry()
        myGeometry.addPoint([0,1,1])
        points = [[0,1,0], [1,1,2], [3,54,1]]
        myGeometry.addPoints(points)
        self.assertTrue(myGeometry.getPoint(0) == [0,1,1])
        self.assertTrue(myGeometry.pointsCount() == 4)


class TestGGeometryLines(unittest.TestCase):

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

    def test_addLineSegment(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addLineSegment([0,1])
        self.assertTrue(G.lineSegmentsCount() == 1)

    def test_addLineSegment_missingPoint(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        try:
            G.addLineSegment([0,10])
        except AssertionError:
            self.assertTrue(True)
            return

        self.assertTrue(False)

    def test_remove_LineSegment(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addLineSegment([0,1])
        G.removeLineSegment(0)
        self.assertTrue(G.lineSegmentsCount() == 0)
        try:
            G.removeLineSegment(0)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_getLineSegment(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addLineSegment([0,1])
        self.assertTrue(G.lineSegmentsCount() == 1)
        line = G.getLineSegment(0)
        self.assertTrue(line == [0,1])


    def test_get_nonexistant_LineSegment(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addLineSegment([0,1])
        self.assertTrue(G.lineSegmentsCount() == 1)
        try:
            G.getLineSegment(1)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_add_lines_dictionary(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        lines = {0: [0,1], 1:[1,2], 2:[2,3]}
        G.addLineSegmentsFromDictionary(lines)
        self.assertTrue(G.lineSegmentsCount() == 3)



class TestGGeometryBSplines(unittest.TestCase):

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

    def test_addBspline(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addBspline([0,1,2])
        self.assertTrue(G.bsplinesCount() == 1)

    def test_addBspline_missingPoint(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        try:
            G.addBspline([0,10,1])
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_remove_Bspline(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addBspline([0,1,2])
        G.removeBspline(0)
        self.assertTrue(G.bsplinesCount() == 0)
        try:
            G.removeLineSegment(0)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_getBspline(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addBspline([0,1,2])
        self.assertTrue(G.bsplinesCount() == 1)
        line = G.getBspline(0)
        self.assertTrue(line == [0,1,2])


    def test_get_nonexistant_Bspline(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        G.addBspline([0,1,2])
        self.assertTrue(G.bsplinesCount() == 1)
        try:
            G.getBspline(1)
        except AssertionError:
            self.assertTrue(True)
            return
        self.assertTrue(False)

    def test_add_Bspline_dictionary_nonexisting(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        lines = {0: [0,1,4], 1:[1,2,3]}
        try:
            G.addBsplinesFromDictionary(lines)
        except AssertionError:
            self.assertTrue(True)
        self.assertTrue(False)

    def test_add_Bspline_dictionary_nonexisting(self):
        #Try importing qmesh. If an ImportError exception is thrown, then
        # the test might be run in a development setting, and see if qmesh
        # can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.append(self.project_root_path)
            import qmesh3
        qmesh3.LOG.setLevel('WARNING')
        G = qmesh3.mesh.Geometry()
        G.addPoints([[0,0,0], [0,1,0], [1,1,0], [1,0,0]])
        lines = {0: [0,1,2], 1:[1,2,3]}
        try:
            G.addBsplinesFromDictionary(lines)
        except AssertionError:
            self.assertTrue(False)
        self.assertTrue(G.bsplinesCount() == 2)


if __name__ == '__main__':
    unittest.main()
