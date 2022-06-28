#    Copyright (C) 2016 Alexandros Avdis and others. See the AUTHORS file for a full list of copyright holders.
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

class Branch(object):
    '''Class Branch, associating a root-node to barches emanating from that node.'''
    def __init__(self, rootNode=None, subBranches=[]):
        '''Class initialisation method, returns an "empty" branch by default, where
        the root node is "None", associated to an empty list of sub-braches.
        :arg rootNode: The node the branch emanates from
        :arg list subBranches: a list of branch-type objects, emanating from the branch end-node
        '''
        self.rootNode = rootNode
        self.subBranches = []#subBranches

class Tree(object):
    '''Class Tree, representing connected branches.'''
    def __init__(self, rootBranch=Branch()):
        '''Class initialisation method, returns an "empty" tree by default, whithout
        any branches. Given the way braches are encoded, only a reference to the "root"
        branch need be stored.
        :arg Branch rootBranch: The root branch of the tree.
        '''
        self.rootBranch = rootBranch
    def nodes(self):
        '''Iterator through nodes of the tree.'''
        branches = [self.rootBranch]
        while len(branches) > 0:
            subBranches = []
            for branch in branches:
                subBranches.extend(branch.subBranches)
                yield branch.rootNode
            branches = subBranches
    def nodesWithChildren(self):
        '''Iterator though nodes of the tree, method also yields the children nodes.'''
        branches = [self.rootBranch]
        while len(branches) > 0:
            subBranches = []
            for branch in branches:
                children = []
                for subBranch in branch.subBranches:
                    children.append(subBranch.rootNode)
                    subBranches.append(subBranch)
                yield (branch.rootNode, children)
            branches = subBranches
    def branchesWithParents(self):
        '''Iterator through branches of the tree, method also yields the parent branch.'''
        yield (self.rootBranch, None)
        branches = [self.rootBranch]
        while len(branches) > 0:
            parentBranches = []
            childBranches = []
            for branch in branches:
                for subBranch in branch.subBranches:
                    parentBranches.append(branch)
                    childBranches.append(subBranch)
            for branch,parent in zip(childBranches, parentBranches):
                yield (branch, parent)
            branches = childBranches
    def nodeList(self):
        '''Returns list of all tree nodes.
        :returns: List of all tree nodes.
        '''
        nodeList = []
        for node in self.nodes():
            nodeList.append(node)
        return nodeList
    def hasNode(self, node):
        '''Method for checking if tree has given node
        :rtype: bool
        :returns: True if node was found in tree, False otherwise
        '''
        if node in self.nodeList():
            return True
        else:
            return False
    def addNode(self, newNode, parentNode):
        '''Method for adding a new node to the tree, as a child to a given parent node.
        :arg newNode: New node to add.
        :parentNode: The parent node.
        '''
        if parentNode not in self.nodeList():
            errorMsg = 'Could not insert node, the parent node was not found in the tree.'
            raise Exception(errorMsg)
        for (branch, parentBranch) in self.branchesWithParents():
            if parentNode == branch.rootNode:
                newBranch = Branch(rootNode = newNode)
                branch.subBranches.append(newBranch)
                break

    def toSubTreeList(self):
        '''This Method returns a list containing [sub-tree, parent node]. Each sub-tree has one and only one node of the original tree as its root-node, and all branches emanating from it. In this way this method also allows to iterate trhough the nodes of a tree, whilist also supplying parent and child information.'''
        subTreeList = [(self, self.rootBranch.rootNode)]
        subBranches = self.rootBranch.subBranches
        while len(subBranches)>0:
            subsubBranches = []
            for subBranch in subBranches:
                subsubBranches.extend(subBranch.subBranches)
                subTreeList.append([Tree(rootBranch=subBranch), subBranch.rootNode])
            subBranches = subsubBranches
        return subTreeList
    def prune(self, node):
        '''Method for deleting sub-trees from a tree.
        :arg node: Node to prune the tree at
        :rtype: Tree
        :returns: The pruned sub-tree as a separate tree
        '''
        import copy
        if node == self.rootBranch.rootNode:
            removedSubTree = copy.copy(self)
            self.rootBranch=Branch()
            return removedSubTree
        for (branch,parentBranch) in self.branchesWithParents():
            if branch.rootNode == node:
                removedSubTree = Tree(branch)
                parentBranch.subBranches.remove(branch)
                return removedSubTree
    def getRootNode(self):
        '''Returns the root node of the tree
        :rtype: Tree
        :returns: Tree root node
        '''
        return self.rootBranch.rootNode
    def isEmpty(self):
        '''Check it current tree is empty
        :rtype: bool
        :returns: True if tree is empty, false if not.
        '''
        if self.getRootNode() == None:
            return True
        else:
            return False


class GridQuadTree(object):
    def __init__(self, xiData, etaData):
        #construct a list of lists, containing the xi-indices that demarcate the
        #regions at each level of the quad tree.
        xiLists = []
        step = 1
        while True:
            xiIndices = list(range(0, len(xiData[0,:]), step))
            if xiIndices[-1] != len(xiData[0,:])-1:
                xiIndices.append(len(xiData[0,:])-1)
            step *= 2
            xiLists.append(xiIndices)
            if len(xiIndices)==2:
                break
        #Construct the xi-tree
        self.xiTree = Tree(Branch(rootNode=xiLists.pop()))
        parentNodes = [self.xiTree.getRootNode()]
        while len(xiLists)!=0:
            childNodes = []
            parentNodesCounter=0
            parentNode = parentNodes[parentNodesCounter]
            xiIndices = xiLists.pop()
            for xiCounter, xiIndex in enumerate(xiIndices):
                if xiIndex == len(xiData[0,:])-1:
                    parentNodes = childNodes
                    break
                if xiIndices[xiCounter+1] > parentNode[1]:
                    parentNodesCounter += 1
                    parentNode = parentNodes[parentNodesCounter]
                childNode = [xiIndex, xiIndices[xiCounter+1]]
                childNodes.append(childNode)
                if childNode == parentNode:
                    pass
                else:
                    self.xiTree.addNode(childNode, parentNode)
        #construct a list of lists, containing the eta-indices that demarcate the
        #regions at each level of the quad tree.
        etaLists = []
        step = 1
        while True:
            etaIndices = list(range(0, len(etaData[:,0]), step))
            if etaIndices[-1] != len(etaData[:,0])-1:
                etaIndices.append(len(etaData[:,0])-1)
            step *= 2
            etaLists.append(etaIndices)
            if len(etaIndices)==2:
                break
        #construct the eta-tree
        self.etaTree = Tree(Branch(rootNode=etaLists.pop()))
        parentNodes = [self.etaTree.getRootNode()]
        while len(etaLists) !=0:
            childNodes = []
            parentNodesCounter=0
            parentNode = parentNodes[parentNodesCounter]
            etaIndices = etaLists.pop()
            for etaCounter, etaIndex in enumerate(etaIndices):
                if etaIndex == len(etaData[:,0])-1:
                    parentNodes = childNodes
                    break
                if etaIndices[etaCounter+1] > parentNode[1]:
                    parentNodesCounter += 1
                    parentNode = parentNodes[parentNodesCounter]
                childNode = [etaIndex, etaIndices[etaCounter+1]]
                childNodes.append(childNode)
                if childNode == parentNode:
                    pass
                else:
                    self.etaTree.addNode(childNode, parentNode)
        #Store references to point coordinates data.
        self.xiData = xiData
        self.etaData = etaData

    def locate_xiCoord(self, xiCoord):
        '''Locate the abscissa coords stored in the tree surrounding the input coordinate.
        :arg xiCoord: coordinate along abscissa to locate its surrounding grid points
        :rtype: float or None
        :returns: List containing the indices of the surrounding points on the abscissa.
                  If the given coordinate exist in the tree grid, the index of the grid point is returned rather than a pair of indices.
                  None is rerurned for points outside the tree coverage
        '''
        branch = self.xiTree.rootBranch
        node = branch.rootNode
        #If point is outside tree, return None
        if not(self.xiData[0,node[0]] <= xiCoord <= self.xiData[0,node[1]]):
            return None
        else:
            while len(branch.subBranches) != 0:
                for subBranch in branch.subBranches:
                    node = subBranch.rootNode
                    if self.xiData[0,node[0]] <= xiCoord <= self.xiData[0,node[1]]:
                        branch = subBranch
                        break
            if self.xiData[0,node[0]] == xiCoord:
                return [node[0]]
            elif self.xiData[0,node[1]] == xiCoord:
                return [node[1]]
            else:
                return node

    def locate_etaCoord(self, etaCoord):
        '''Locate the ordinate coords stored in the tree surrounding the input coordinate.
        :arg etaCoord: coordinate along ordinate to locate its surrounding grid points
        :rtype: float or None
        :returns: List containing the indices of the surrounding points on the ordinate.
                  If the given coordinate exist in the tree grid, the index of the grid point is returned rather than a pair of indices.
                  None is rerurned for points outside the tree coverage
        '''
        branch = self.etaTree.rootBranch
        node = branch.rootNode
        #If point is outside tree, return None
        if not(self.etaData[node[0],0] <= etaCoord <= self.etaData[node[1],0]):
            return None
        else:
            while len(branch.subBranches) != 0:
                for subBranch in branch.subBranches:
                    node = subBranch.rootNode
                    if self.etaData[node[0],0] <= etaCoord <= self.etaData[node[1],0]:
                        branch = subBranch
                        break
            if self.etaData[node[0],0] == etaCoord:
                return [node[0]]
            elif self.etaData[node[1],0] == etaCoord:
                return [node[1]]
            else:
                return node

    def locatePoint(self, point):
        '''Locate the coordinates stored in the tree surrounding the given point.
        :arg point: Coordinates list of the point to be located.
        :rtype: float or None
        :returns: List containing two lists: the first list contains the indices of the surrounding points along the abscissa
                  and the second contains the indices of the surrounding points along the ordinate.
                  If the given coordinate exist in the tree grid, the index of the grid point is returned rather than a pair of indices.
                  None is rerurned for points outside the tree coverage
        '''
        xiCoord = point[0]
        etaCoord = point[1]
        surrounding_xiIndices = self.locate_xiCoord(xiCoord)
        surrounding_etaIndices = self.locate_etaCoord(etaCoord)
        if surrounding_xiIndices != None and surrounding_etaIndices!= None:
            return [surrounding_xiIndices, surrounding_etaIndices] 
        else:
            return None
