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
# so we import local python before any other
import sys
sys.path.insert(0,"../../../")
sys.path.insert(0,"../")
import qmesh

class TestGlobe(unittest.TestCase):
    '''Test qmesh functionality on idealised shapes in planet-centered-cartesian.'''

    def setUp(self):
            qmesh.LOG.setLevel('ERROR')
            import os
            self.thisPath = os.path.dirname(os.path.realpath(__file__))
            self.equatorLine_filename = self.thisPath+'/equator.shp'
            self.islandShoreline_filename = self.thisPath+'/circularIsland.shp'

    def create_equator_shapefiles(self):
        '''Test primitive shape creation utilities. Create a line representative of the equator, in EPSG:4326, and write to a shapefile.'''
        #Try equator in line-shapefiles
        try:
            equatorLine = qmesh.vector.loxodromicLine( \
                                startPoint = (-179.,0.), \
                                trueNorthBearing = 90., \
                                numbPoints=20, \
                                loopArounds = 0, \
                                endPoint = (179.,0.), \
                                coordRefSystem_string = "EPSG:4326")
            equatorLine.asShapes().writeFile(self.equatorLine_filename)
        except AssertionError:
            self.assertTrue(False)

    def create_circular_island(self):
        '''Test primitive shape creation utilities. Create a circle, in EPSG:4326, and write to a shapefile.'''
        try:
            islandShoreline = qmesh.vector.Circle(
                                centerPointXi = 0.0,
                                centerPointEta = -20.0,
                                radius = 10,
                                numbPoints = 20,
                                coordRefSystem_string = "EPSG:4326")
            islandShoreline.asShapes().writeFile(self.islandShoreline_filename)
        except AssertionError:
            self.assertTrue(False)

    def test_halfGlobe(self):
        '''Test meshing of idealised shapes in planet-centered-cartesian. Meshing the northen hemisphere.'''
        function_name = self.test_halfGlobe.__name__
        #Try initialising qgis API
        try:
            qmesh.initialise()
        except AssertionError:
            self.assertTrue(False)
        #Create shapefiles.
        self.create_equator_shapefiles()
        #read-in Shapefiles
        try:
            boundary = qmesh.vector.Shapes()
            boundary.fromFile(self.equatorLine_filename)
            loopShapes = qmesh.vector.identifyLoops(boundary,
                      isGlobal=True, defaultPhysID=1, fixOpenLoops=True, extraPointsPerVertex=3)
            polygonShapes = qmesh.vector.identifyPolygons(loopShapes,
                      isGlobal=True, meshedAreaPhysID=10)
            raster = qmesh.raster.gradationToShapes()
            raster.setShapes(loopShapes)
            raster.setRasterBounds(-180.0, 180.0, -90.0, 90.0)
            raster.setRasterResolution(1000,500)
            # mesh size in m, distance in deg
            raster.setGradationParameters(10000.0,1000000.0,30.0,1.0)
            raster.calculateLinearGradation()
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loopShapes, polygonShapes)
            domain.setMeshMetricField(raster)
            domain.setTargetCoordRefSystem('PCC', fldFillValue=100000.0)
        except AssertionError:
            self.assertTrue(False)
        #Try meshing. Note, qmesh.gmshDomain() defaults to del2d meshing
        # algo, but in this case it does not give good results.
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo', 
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_islandInGlobe(self):
        '''Test meshing the globe, with circular island.'''
        function_name = self.test_islandInGlobe.__name__
        #Initialise qgis API
        qmesh.initialise()
        #Create shapefiles.
        self.create_circular_island()
        #read-in Shapefiles and set-up domain
        boundary = qmesh.vector.Shapes()
        boundary.fromFile(self.islandShoreline_filename)
        loopShapes = qmesh.vector.identifyLoops(boundary,
                  isGlobal=True, defaultPhysID=1, fixOpenLoops=True)
        polygonShapes = qmesh.vector.identifyPolygons(loopShapes,
                  isGlobal=True, southPoleCoordinates=(0.0,-20.0), meshedAreaPhysID=10)
        raster = qmesh.raster.gradationToShapes()
        raster.setShapes(loopShapes)
        raster.setRasterBounds(-180.0, 180.0, -90.0, 90.0)
        raster.setRasterResolution(800,400)
        raster.setGradationParameters(100000.0,1000000.0,30.0)
        raster.calculateLinearGradation()
        raster.writeNetCDF(self.thisPath+'/'+function_name+'.nc')
        domain = qmesh.mesh.Domain()
        domain.setGeometry(loopShapes, polygonShapes)
        domain.setMeshMetricField(raster)
        domain.setTargetCoordRefSystem('PCC', fldFillValue=1000.0)
        #Generate mesh
        domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo', 
                    fldFilename=self.thisPath+'/'+function_name+'.fld',
                    mshFilename=self.thisPath+'/'+function_name+'.msh',
                    gmshAlgo = 'meshadapt')

    def test_loxodromeBoundRegion(self):
        '''Test meshing of idealised shapes in planet-centered-cartesian. Meshing the region bound between two loxodromes.'''
        import numpy as np
        #Try initialising qgis API
        try:
            qmesh.initialise()
        except AssertionError:
            self.assertTrue(False)
        #Create shapefiles.
        ldrome_1 = qmesh.vector.loxodromicLine((0.0, -85.0), 50.0,100, 0, (np.nan, 85.0), "EPSG:4326")
        ldrome_2 = qmesh.vector.loxodromicLine((0.0, -85.0), 40.0,100, 0, (np.nan, 85.0), "EPSG:4326")
        ldrome_1.asShapes().writeFile(self.thisPath+'/loxodrome1.shp')
        ldrome_2.asShapes().writeFile(self.thisPath+'/loxodrome2.shp')

suite = unittest.TestLoader().loadTestsFromTestCase(TestGlobe)


if __name__ == '__main__':
    unittest.main()
