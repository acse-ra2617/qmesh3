#    Copyright (C) 2016 Alexandros Avdis and others.
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

#class TestBranch(unittest.TestCase):
#    def testEmptyBranch(self):
#        emptyBranch = qmesh.lib.Trees.Branch()
#        self.assertIsNone(emptyBranch.rootNode)
#        self.assertEqual(emptyBranch.subBranches, [])
#
#class TestTree(unittest.TestCase):
#    def testEmptyTree(self):
#        emptyBranch = qmesh.lib.Trees.Branch()
#        emptyTree = qmesh.lib.Trees.Tree(emptyBranch)

class TestGridQuadTree(unittest.TestCase):
    '''Testing the GridQuadTree class in the lib module of qmesh.'''
    def setUp(self):
        '''Method setting-up test fixture.

        This method is not aimed to be invoked by the user. Various resources needed and
        produced by the test are defined here:
        1) The path of the testing function and the project root path (qmesh).
        2) Basic data-structures used in testing.
        3) The correct results, we expect the test-cases to return.
        '''
        import os
        import sys
        import numpy
        #Construct absolute path to this test file, without the test file-name
        self.this_path = os.path.dirname(os.path.realpath(__file__))
        #Construct absolute path to the test-type directory
        self.project_test_type_path = os.path.split(self.this_path)[0]
        #Construct absolute path to the tests directory
        self.project_tests_path = os.path.split(self.project_test_type_path)[0]
        #Construct absolute path to the project root
        self.project_root_path = os.path.split(self.project_tests_path)[0]
        #Try importing qmesh. If an ImportError exception is thrown, then the test might be run
        # in a development setting, and see if qmesh can be found at the project root.
        try:
            import qmesh3
        except ImportError:
            sys.path.insert(0,self.project_root_path)
            import qmesh3
        #Define a mesh-grid, by calculating the point coordinates, and populate a qmesh
        # GridQuadTree object
        xi = numpy.linspace(0,10,11)
        eta = numpy.linspace(-10,10,21)
        xiData, etaData = numpy.meshgrid(xi, eta)
        self.gridQuadTree = qmesh3.lib.Trees.GridQuadTree(xiData, etaData)
        #Define the expected resutls of the quad tree algorithm on the test mesh-grid.
        self.expected_xiNodes = [[0,10],
                                 [0,8],[8,10],
                                 [0,4],[4,8],[8,9],[9,10],
                                 [0,2],[2,4],[4,6],[6,8],
                                 [0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]
        self.expected_etaNodes = [[0,20],
                                  [0,16],[16,20],
                                  [0,8],[8,16],[16,18],[18,20],
                                  [0,4],[4,8],[8,12],[12,16],[16,17],[17,18],[18,19],[19,20],
                                  [0,2],[2,4],[4,6],[6,8],[8,10],[10,12],[12,14],[14,16],
                                  [0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],
                                               [10,11],[11,12],[12,13],[13,14],[14,15],[15,16]]
    def testExpectedNodes(self):
        '''Testing the quad-tree node calculation.

        This method tests if the calculation of the array indices where partitions are drawn in a
        structured quad-grid are correct.
        '''
        for xiNode, expected_xiNode in zip(self.gridQuadTree.xiTree.nodes(), self.expected_xiNodes):
            self.assertEqual(xiNode, expected_xiNode)
        for etaNode, expected_etaNode in zip(self.gridQuadTree.etaTree.nodes(), self.expected_etaNodes):
            self.assertEqual(etaNode, expected_etaNode)
    def test_locate_xiCoord(self):
        self.assertIsNone(self.gridQuadTree.locate_xiCoord(-0.1))
        self.assertEqual(self.gridQuadTree.locate_xiCoord(0.0), [0])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(0.1), [0,1])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(1.1), [1,2])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(2.1), [2,3])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(3.1), [3,4])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(4.1), [4,5])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(5.1), [5,6])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(6.1), [6,7])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(7.1), [7,8])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(8.1), [8,9])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(9.1), [9,10])
        self.assertEqual(self.gridQuadTree.locate_xiCoord(10.0), [10])
        self.assertIsNone(self.gridQuadTree.locate_xiCoord(10.1))
    def test_locate_etaCoord(self):
        self.assertIsNone(self.gridQuadTree.locate_etaCoord(-10.1))
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-10.), [0])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-9.1), [0,1])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-8.1), [1,2])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-7.1), [2,3])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-6.1), [3,4])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-5.1), [4,5])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-4.1), [5,6])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-3.1), [6,7])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-2.1), [7,8])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-1.1), [8,9])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(-0.1), [9,10])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(0.0), [10])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(0.9), [10,11])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(1.9), [11,12])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(2.9), [12,13])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(3.9), [13,14])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(4.9), [14,15])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(5.9), [15,16])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(6.9), [16,17])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(7.9), [17,18])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(8.9), [18,19])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(9.9), [19,20])
        self.assertEqual(self.gridQuadTree.locate_etaCoord(10), [20])
        self.assertIsNone(self.gridQuadTree.locate_etaCoord(10.9))
    def test_locatePoint(self):
        self.assertIsNone(self.gridQuadTree.locatePoint([-0.1,-10.1]))
        self.assertEqual(self.gridQuadTree.locatePoint([0.0,-10.]), [[0],[0]])
        self.assertEqual(self.gridQuadTree.locatePoint([0.1,-10.]), [[0,1],[0]])
        self.assertEqual(self.gridQuadTree.locatePoint([0.,-9.1]), [[0],[0,1]])
        self.assertEqual(self.gridQuadTree.locatePoint([0.1,-9.1]), [[0,1],[0,1]])
        self.assertEqual(self.gridQuadTree.locatePoint([1.1,-8.1]), [[1,2],[1,2]])
        self.assertEqual(self.gridQuadTree.locatePoint([2.1,-7.1]), [[2,3],[2,3]])
        self.assertEqual(self.gridQuadTree.locatePoint([3.1,-6.1]), [[3,4],[3,4]])
        self.assertEqual(self.gridQuadTree.locatePoint([4.1,-5.1]), [[4,5],[4,5]])
        self.assertEqual(self.gridQuadTree.locatePoint([5.1,-4.1]), [[5,6],[5,6]])
        self.assertEqual(self.gridQuadTree.locatePoint([6.1,-3.1]), [[6,7],[6,7]])
        self.assertEqual(self.gridQuadTree.locatePoint([7.1,-2.1]), [[7,8],[7,8]])
        self.assertEqual(self.gridQuadTree.locatePoint([8.1,-1.1]), [[8,9],[8,9]])
        self.assertEqual(self.gridQuadTree.locatePoint([9.1,-0.1]), [[9,10],[9,10]])
        self.assertEqual(self.gridQuadTree.locatePoint([10.,0]), [[10],[10]])
        self.assertEqual(self.gridQuadTree.locatePoint([9.1,0]), [[9,10],[10]])
        self.assertIsNone(self.gridQuadTree.locatePoint([11.1,0]))
        self.assertEqual(self.gridQuadTree.locatePoint([10,-0.1]), [[10],[9,10]])
        self.assertEqual(self.gridQuadTree.locatePoint([10,0.9]), [[10],[10,11]])
        self.assertEqual(self.gridQuadTree.locatePoint([9.1,0.9]), [[9,10],[10,11]])
        self.assertEqual(self.gridQuadTree.locatePoint([8.1,1.9]), [[8,9],[11,12]])
        self.assertEqual(self.gridQuadTree.locatePoint([7.1,2.9]), [[7,8],[12,13]])
        self.assertEqual(self.gridQuadTree.locatePoint([6.1,3.9]), [[6,7],[13,14]])
        self.assertEqual(self.gridQuadTree.locatePoint([5.1,4.9]), [[5,6],[14,15]])
        self.assertEqual(self.gridQuadTree.locatePoint([4.1,5.9]), [[4,5],[15,16]])
        self.assertEqual(self.gridQuadTree.locatePoint([3.1,6.9]), [[3,4],[16,17]])
        self.assertEqual(self.gridQuadTree.locatePoint([2.1,7.9]), [[2,3],[17,18]])
        self.assertEqual(self.gridQuadTree.locatePoint([1.1,8.9]), [[1,2],[18,19]])
        self.assertEqual(self.gridQuadTree.locatePoint([0.1,9.9]), [[0,1],[19,20]])
        self.assertEqual(self.gridQuadTree.locatePoint([0.,10.]), [[0],[20]])
        self.assertIsNone(self.gridQuadTree.locatePoint([-0.9,10.]))

suite = unittest.TestLoader().loadTestsFromTestCase(TestGridQuadTree)

if __name__ == '__main__':
    unittest.main()
