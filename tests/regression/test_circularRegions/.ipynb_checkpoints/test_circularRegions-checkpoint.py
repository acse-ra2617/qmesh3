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
import numpy as np
import qgis.core
# so we import local python before any other
import sys
# if running the test manually
sys.path.insert(0,"../../../")
# if running the test via Makefile
sys.path.insert(0,"../")
import qmesh3

class TestCircularRegions(unittest.TestCase):
    '''Test region insertion in qmesh.'''

    def setUp(self):
        qmesh.LOG.setLevel('WARNING')
        import os
        self.thisPath = os.path.dirname(os.path.realpath(__file__))
        self.circularDomain_line_filename = self.thisPath+'/test_outerCircle_lines.shp'
        self.circularDomain_polygon_filename = self.thisPath+'/test_outerCircle_polygons.shp'
        self.circularRegion30_line_filename = self.thisPath+'/test_circularRegion30_lines.shp'
        self.circularRegion30_polygon_filename = \
                                          self.thisPath+'/test_circularRegion30_polygons.shp'
        self.circularRegions5_line_filename = \
                                  self.thisPath+'/test_circularRegions_5degRadius_lines.shp'
        self.circularRegions5_polygon_filename = \
                              self.thisPath+'/test_circularRegions_5degRadius_polygons.shp'
        self.create_outerRegion_shapefiles()
        self.create_circularRegion30_shapefiles()
        self.create_circularRegions5_shapefiles()


    def create_outerRegion_shapefiles(self):
        '''Test shape creation with qmesh API. Create a circle, in EPSG:4326, centered at (0,0) and radius of 80. Write to files as line and polygon.'''
        #Try creating line feature.
        try:
            #Calculate point coordinates
            numberPoints=1000
            deltaAngle=360./numberPoints
            centerPointCoords = [0.,0.]
            circleRadius = 80.0
            outer_boundary = qmesh.vector.Circle(centerPointCoords[0], centerPointCoords[1],
                    circleRadius, numberPoints, "EPSG:4326")
        except AssertionError:
            self.assertTrue(False)
        # Try creating a polygon from the line
        try:
            # Create line feature
            loops = qmesh.vector.identifyLoops(outer_boundary.asShapes(),
                    defaultPhysID=1000, fixOpenLoops=True)
        except AssertionError:
            self.assertTrue(False)
        # Write shapefiles
        try:
            loops.writeFile(self.circularDomain_line_filename)
        except AssertionError:
            self.assertTrue(False)

    def create_circularRegion30_shapefiles(self):
        '''Test shape creation with qmesh API. Create a circle, in EPSG:4326, centered at (0,10) and radius of 30. Write to files as line and polygon.'''
        #Try greating line feature.
        try:
            #Calculate point coordinates
            numberPoints=1000
            deltaAngle=360./numberPoints
            centerPointCoords = [0.,10.]
            circleRadius = 30.0
            boundary = qmesh.vector.Circle(centerPointCoords[0], centerPointCoords[1],
                    circleRadius, numberPoints, "EPSG:4326")
        except AssertionError:
            self.assertTrue(False)
        #Try creating a polygon from the line
        try:
            # Create line feature
            loops = qmesh.vector.identifyLoops(boundary.asShapes(),
                    defaultPhysID=1000, fixOpenLoops=True)
        except AssertionError:
            self.assertTrue(False)
        # Write shapefiles
        try:
            loops.writeFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)


    def create_circularRegions5_shapefiles(self):
        '''Test shape creation with qmesh API. Create a circle, in EPSG:4326, centered at (20,0) and radius of 5. Write to files as line and polygon.'''
        #Try greating line features.
        try:
            circles = qmesh.vector.Shapes()
            circles.setCoordRefSystemFromString('EPSG:4326')
            circles.setShapeType(qgis.core.QgsWkbTypes.LineString)
            #
            numberPoints=1000
            deltaAngle=360./numberPoints
            circleRadius = 5.0
            # Circle 1
            centerPointCoords = [20.,0.]
            circle1 = qmesh.vector.Circle(centerPointCoords[0], centerPointCoords[1],
                    circleRadius, numberPoints, "EPSG:4326")
            circles.addFeature(circle1.asQgsFeature())
            # Circle 2
            centerPointCoords = [0.,10.]
            circle2 = qmesh.vector.Circle(centerPointCoords[0], centerPointCoords[1],
                    circleRadius, numberPoints, "EPSG:4326")
            circles.addFeature(circle2.asQgsFeature())
            # Circle 3
            circles_loops = qmesh.vector.identifyLoops(circles,
                      defaultPhysID=1000, fixOpenLoops=True)
            circles_polygons = qmesh.vector.identifyPolygons(circles_loops,
                      meshedAreaPhysID = 80000)
        except AssertionError:
            self.assertTrue(False)
        # Write shapefiles
        try:
            circles_loops.writeFile(self.circularRegions5_line_filename)
        except AssertionError:
            self.assertTrue(False)


    def test_singleRegion_PCC_noMeshMetric(self):
        '''Test reading-in lines and polygons from file and generating a mesh, in planet-centered-cartesian.

        An incremental step towards testing
        "region insertion" functionality in qmesh. However, no regions are
        created nor inserted in this test. Only the "outer" domain is created
        and meshed. The input domain is in EPSG:4326 and the output mesh is
        in 'Planet-Centered-Cartesian'. The domain is a circle centered at
        (0.0, 0.0)degrees, with an 80 degree radius.'''
        function_name = self.test_singleRegion_PCC_noMeshMetric.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct lines and polygons
        try:
            #Create polygon feature
            loops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
            polygons = qmesh.vector.identifyPolygons(loops,
                    meshedAreaPhysID = 10000)
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setTargetCoordRefSystem('PCC')
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_singleRegion_EPSG6933_noMeshMetric_deln(self):
        '''Test reading-in lines and polygons from file and generating a mesh, in EPSG:6933.

        Rationale: This function is an incremental step towards testing
        "region insertion" functionality in qmesh. However, no regions are
        created nor inserted in this test, only the "outer" domain is created
        and meshed. The input domain is in EPSG:4326 and the output mesh is
        in EPSG:6933. The domain is a circle centered at (0.0, 0.0)degrees,
        with an 80 degree radius. The gmsh Delaunay mesh generation algorithm
        is used'''
        function_name = self.test_singleRegion_EPSG6933_noMeshMetric_deln.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct lines and polygons
        try:
            #Create polygon feature
            loops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
            polygons = qmesh.vector.identifyPolygons(loops,
                    meshedAreaPhysID = 10000)
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setTargetCoordRefSystem('EPSG:6933')
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'del2d')
        except AssertionError:
            self.assertTrue(False)

    def test_singleRegion_EPSG6933_noMeshMetric_frnt(self):
        '''Test reading-in lines and polygons from file and generating a mesh, in EPSG:6933.

        Rationale: This function is an incremental step towards testing 
        "region insertion" functionality in qmesh. However no regions are
        created nor inserted in this test, only the "outer" domain is created
        and meshed. The input domain is in EPSG:4326 and the output mesh is
        in EPSG:6933. The  domain is a circle centered at (0.0, 0.0)degrees,
        with an 80 degree radius. The gmsh frontal mesh generation algorithm
        is used.'''
        function_name = self.test_singleRegion_EPSG6933_noMeshMetric_frnt.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct lines and polygons
        try:
            #Create polygon feature
            loops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
            polygons = qmesh.vector.identifyPolygons(loops,
                    meshedAreaPhysID = 10000)
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setTargetCoordRefSystem('EPSG:6933')
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_twoRegion_PCC_noMeshMetric(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of two regions, and a mesh is generated in
        Planet-Centered-Cartesian. No mesh gradation is specified.'''
        function_name = self.test_twoRegion_PCC_noMeshMetric.__name__
        # Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        # Insert smaller region into larger.
        try:
            loops, polygons = \
              qmesh.vector.insertRegions(
                outerRegionLoops,\
                outerRegionPolygons,\
                circularRegion30Loops,\
                circularRegion30Polygons)
        except AssertionError:
            self.assertTrue(False)
        # Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setTargetCoordRefSystem('PCC')
        except AssertionError:
            self.assertTrue(False)
        # Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_twoRegion_EPSG6933_noMeshMetric(self):
        '''Test multiple region meshing and labelling.
        
        Domain is composed of two regions, and a mesh is generated 
        in EPSG:6933. No mesh gradation is specified.'''

        function_name = self.test_twoRegion_EPSG6933_noMeshMetric.__name__
        # Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        # Insert smaller region into larger.
        try:
            loops, polygons = \
              qmesh.vector.insertRegions(
                outerRegionLoops,\
                outerRegionPolygons,\
                circularRegion30Loops,\
                circularRegion30Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setTargetCoordRefSystem('EPSG:6933')
        except AssertionError:
            self.assertTrue(False)
        # Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_twoRegion_PCC_gradatedMesh(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of two regions, and a mesh is generated
        in Planet-Centered-Cartesian. Mesh gradates towards inner
        region.'''
        function_name = self.test_twoRegion_PCC_gradatedMesh.__name__
        # Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        # Insert smaller region into larger.
        try:
            loops, polygons = \
                qmesh.vector.insertRegions(
                   outerRegionLoops, outerRegionPolygons,\
                   circularRegion30Loops, circularRegion30Polygons)
        except AssertionError:
            self.assertTrue(False)
        # Create gradated mesh metric raster for inserted region.
        try:
            raster = qmesh.raster.gradationToShapes()
            raster.setShapes(circularRegion30Polygons)
            raster.setRasterBounds(-85.0,85.0,-85.0,85.0)
            raster.setRasterResolution(100,100)
            raster.setGradationParameters(100000.0,500000.0,20.0)
            raster.calculateLinearGradation()
        except AssertionError:
            self.assertTrue(False)
        # Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(raster)
            domain.setTargetCoordRefSystem('PCC')
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'del2d')
        except AssertionError:
            self.assertTrue(False)

    def test_twoRegion_EPSG6933_gradatedMesh_deln(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of two regions, and a mesh is generated
        in EPSG:6933. Mesh gradates towards inner region. The gmsh
        Delaunay mesh generation algorith is used.
        '''
        function_name = self.test_twoRegion_EPSG6933_gradatedMesh_deln.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        #Insert smaller region into larger.
        try:
            loops, polygons = \
                qmesh.vector.insertRegions(
                   outerRegionLoops, outerRegionPolygons,\
                   circularRegion30Loops, circularRegion30Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create gradated mesh metric raster for inserted region.
        try:
            raster = qmesh.raster.gradationToShapes()
            raster.setShapes(circularRegion30Polygons)
            raster.setRasterBounds(-85.0,85.0,-85.0,85.0)
            raster.setRasterResolution(100,100)
            raster.setGradationParameters(1.0,50.0,20.0)
            raster.calculateLinearGradation()
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(raster)
            domain.setTargetCoordRefSystem('EPSG:6933', fldFillValue=1000.0)
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'del2d')
        except AssertionError:
            self.assertTrue(False)

    def test_twoRegion_EPSG6933_gradatedMesh_frnt(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of two regions, and a mesh is generated
        in EPSG:6933. Mesh gradates towards inner region.
        The gmsh 2D frontal mesh generation algorithm is used.
        '''
        function_name = self.test_twoRegion_EPSG6933_gradatedMesh_frnt.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        #Insert smaller region into larger.
        try:
            loops, polygons = \
                qmesh.vector.insertRegions(
                   outerRegionLoops, outerRegionPolygons,\
                   circularRegion30Loops, circularRegion30Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create gradated mesh metric raster for inserted region.
        try:
            raster = qmesh.raster.gradationToShapes()
            raster.setShapes(circularRegion30Polygons)
            raster.setRasterBounds(-85.0,85.0,-85.0,85.0)
            raster.setRasterResolution(100,100)
            raster.setGradationParameters(1.0,50.0,20.0)
            raster.calculateLinearGradation()
        except AssertionError:
            self.assertTrue(False)
        #Create Gmsh domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(raster)
            domain.setTargetCoordRefSystem('EPSG:6933', fldFillValue=1000.0)
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)

    def test_multipleRegion_PCC(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of four regions, and a mesh is generated in
        Planet-Centered-Cartesian.
        '''
        function_name = self.test_multipleRegion_PCC.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
            circularRegion5Lines = qmesh.vector.Shapes()
            circularRegion5Lines.fromFile(self.circularRegions5_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        circularRegion5Loops = qmesh.vector.identifyLoops(circularRegion5Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion5Polygons = qmesh.vector.identifyPolygons(circularRegion5Loops,
                    meshedAreaPhysID=50000)
        #Insert smaller region into larger.
        try:
            loops, polygons = \
               qmesh.vector.shapefileTools.insertRegions(
                  outerRegionLoops, outerRegionPolygons,\
                  circularRegion30Loops, circularRegion30Polygons)
            loops, polygons = \
               qmesh.vector.shapefileTools.insertRegions(
                  loops, polygons,\
                  circularRegion5Loops, circularRegion5Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create gradated mesh metric raster for inserted region.
        try:
            rasterI = qmesh.raster.gradationToShapes()
            rasterI.setShapes(circularRegion5Polygons)
            rasterI.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterI.setRasterResolution(200,200)
            rasterI.setGradationParameters(50000.0,500000.0,30.0)
            rasterI.calculateLinearGradation()
            rasterII = qmesh.raster.gradationToShapes()
            rasterII.setShapes(circularRegion30Polygons)
            rasterII.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterII.setRasterResolution(200,200)
            rasterII.setGradationParameters(100000.0,500000.0,20.0)
            rasterII.calculateLinearGradation()
            rasterIII = qmesh.raster.minimumRaster([rasterI, rasterII])
        except AssertionError:
            self.assertTrue(False)
        #Create domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(rasterIII)
            domain.setTargetCoordRefSystem('PCC', fldFillValue=0.0)
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'del2d')
        except AssertionError:
            self.assertTrue(False)

    def test_multipleRegion_EPSG6933_deln(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of four regions, and a mesh is generated in EPSG:6933.
        The gmsh 2D Delaunay mesh generation algorithm is used.
        '''
        function_name = self.test_multipleRegion_EPSG6933_deln.__name__
        #Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
            circularRegion5Lines = qmesh.vector.Shapes()
            circularRegion5Lines.fromFile(self.circularRegions5_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        circularRegion5Loops = qmesh.vector.identifyLoops(circularRegion5Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion5Polygons = qmesh.vector.identifyPolygons(circularRegion5Loops,
                    meshedAreaPhysID=50000)
        #Insert smaller region into larger.
        try:
            loops, polygons = \
               qmesh.vector.insertRegions(
                  outerRegionLoops, outerRegionPolygons,\
                  circularRegion30Loops, circularRegion30Polygons)
            loops, polygons = \
               qmesh.vector.insertRegions(
                  loops, polygons,\
                  circularRegion5Loops, circularRegion5Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create gradated mesh metric raster for inserted region.
        try:
            rasterI = qmesh.raster.gradationToShapes()
            rasterI.setShapes(circularRegion5Polygons)
            rasterI.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterI.setRasterResolution(200,200)
            rasterI.setGradationParameters(2.0,50.0,10.0)
            rasterI.calculateLinearGradation()
            rasterII = qmesh.raster.gradationToShapes()
            rasterII.setShapes(circularRegion30Polygons)
            rasterII.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterII.setRasterResolution(200,200)
            rasterII.setGradationParameters(10.0,50.0,20.0)
            rasterII.calculateLinearGradation()
            rasterIII = qmesh.raster.minimumRaster([rasterI, rasterII])
        except AssertionError:
            self.assertTrue(False)
        #Create domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(rasterIII)
            domain.setTargetCoordRefSystem('EPSG:6933', fldFillValue=0.0)
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'del2d')
        except AssertionError:
            self.assertTrue(False)

    def test_multipleRegion_EPSG6933_frnt(self):
        '''Test multiple region meshing and labelling.

        Domain is composed of four regions, and a mesh is generated in EPSG:6933.
        The Gmsh 2D frontal mesh generation algorithm is used'''
        function_name = self.test_multipleRegion_EPSG6933_frnt.__name__
        # Read-in shapefiles
        try:
            outerRegionLines = qmesh.vector.Shapes()
            outerRegionLines.fromFile(self.circularDomain_line_filename)
            circularRegion30Lines = qmesh.vector.Shapes()
            circularRegion30Lines.fromFile(self.circularRegion30_line_filename)
            circularRegion5Lines = qmesh.vector.Shapes()
            circularRegion5Lines.fromFile(self.circularRegions5_line_filename)
        except AssertionError:
            self.assertTrue(False)
        # Construct loops and polygons
        outerRegionLoops = qmesh.vector.identifyLoops(outerRegionLines,
                    defaultPhysID=1000, fixOpenLoops=True)
        outerRegionPolygons = qmesh.vector.identifyPolygons(outerRegionLoops,
                    meshedAreaPhysID=10000)
        circularRegion30Loops = qmesh.vector.identifyLoops(circularRegion30Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion30Polygons = qmesh.vector.identifyPolygons(circularRegion30Loops,
                    meshedAreaPhysID=10000)
        circularRegion5Loops = qmesh.vector.identifyLoops(circularRegion5Lines,
                    defaultPhysID=1000, fixOpenLoops=True)
        circularRegion5Polygons = qmesh.vector.identifyPolygons(circularRegion5Loops,
                    meshedAreaPhysID=50000)
        # Insert smaller region into larger.
        try:
            loops, polygons = \
               qmesh.vector.insertRegions(
                  outerRegionLoops, outerRegionPolygons,\
                  circularRegion30Loops, circularRegion30Polygons)
            loops, polygons = \
               qmesh.vector.insertRegions(
                  loops, polygons,\
                  circularRegion5Loops, circularRegion5Polygons)
        except AssertionError:
            self.assertTrue(False)
        #Create gradated mesh metric raster for inserted region.
        try:
            rasterI = qmesh.raster.gradationToShapes()
            rasterI.setShapes(circularRegion5Polygons)
            rasterI.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterI.setRasterResolution(200,200)
            rasterI.setGradationParameters(2.0,50.0,10.0)
            rasterI.calculateLinearGradation()
            rasterII = qmesh.raster.gradationToShapes()
            rasterII.setShapes(circularRegion30Polygons)
            rasterII.setRasterBounds(-85.0,85.0,-85.0,85.0)
            rasterII.setRasterResolution(200,200)
            rasterII.setGradationParameters(10.0,50.0,20.0)
            rasterII.calculateLinearGradation()
            rasterIII = qmesh.raster.minimumRaster([rasterI, rasterII])
        except AssertionError:
            self.assertTrue(False)
        #Create domain object
        try:
            domain = qmesh.mesh.Domain()
            domain.setGeometry(loops, polygons)
            domain.setMeshMetricField(rasterIII)
            domain.setTargetCoordRefSystem('EPSG:6933', fldFillValue=0.0)
        except AssertionError:
            self.assertTrue(False)
        #Meshing with Gmsh
        try:
            domain.gmsh(geoFilename=self.thisPath+'/'+function_name+'.geo',
                        fldFilename=self.thisPath+'/'+function_name+'.fld',
                        mshFilename=self.thisPath+'/'+function_name+'.msh',
                        gmshAlgo = 'front2d')
        except AssertionError:
            self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
