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

import qgis.core
import os
import glob
import sys
import numpy as np
import qmesh3
import unittest

def generate_pointsShapefile(points, outputFilename, geoCS_WKT):
    '''Todo: add docstring'''
    import PyQt5.QtCore
    qgis.core.QgsApplication.setPrefixPath('/usr', True)
    qgis.core.QgsApplication.initQgis()

    coordinateReferenceSystem = qgis.core.QgsCoordinateReferenceSystem()
    coordinateReferenceSystem.createFromWkt(geoCS_WKT)
    assert coordinateReferenceSystem.isValid()

    fields =  qgis.core.QgsFields()
    fields.append( qgis.core.QgsField("ID", PyQt5.QtCore.QVariant.Int) )

    writer = qgis.core.QgsVectorFileWriter(outputFilename, "CP1250", fields, qgis.core.QgsWkbTypes.Point, coordinateReferenceSystem, "ESRI Shapefile")
    if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
        print("Error when creating shapefile: ", writer.hasError())

    #
    feature = qgis.core.QgsFeature()
    feature.setFields(fields)
    for point in points:
        QGIS_Point = qgis.core.QgsPoint(point[0],point[1])
        feature.setGeometry(QGIS_Point)
        feature.setAttribute("ID",0)
        writer.addFeature(feature)

    #delete the writer to flush features to disk (optional)
    del writer

class TestProjections(unittest.TestCase):
    '''Todo: add test documentation as class docstring '''

    def setUp(self):
            '''Method setting-up test fixture.'''
            qmesh3.initialise()

    def tearDown(self):
        #Clean-up
        unwantedFiles = glob.glob('test_points.*')
        for unwantedFile in unwantedFiles:
            os.remove(unwantedFile)

    def test_WGS84_2_polarStereographic(self):
        #Create points in shapefile
        sourceCRS_WKT = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'
        generate_pointsShapefile([[0,90],
                                  [-135.0,0],[-90,0],[-45,0],[0,0],[45,0],[90,0],[135.0,0],[180,0],
                                 ],
                                  'test_points.shp', sourceCRS_WKT)
        #Read shapefile
        pointVectorLayer = qmesh3.vector.Shapes()
        pointVectorLayer.fromFile("test_points.shp")
        sourceCRS = pointVectorLayer.getCoordRefSystem()
        self.assertTrue(sourceCRS.isValid())
        #Initialise target projection
        targetCRS = qgis.core.QgsCoordinateReferenceSystem()
        targetCRS.createFromWkt('PROJCS["unnamed",GEOGCS["Normal Sphere (r=0.5)",DATUM["unknown",SPHEROID["sphere",0.5,0]],PRIMEM["Greenwich",0],UNIT["degree",0.01745329251994328]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",90],PARAMETER["central_meridian",90],UNIT["Meter",1]]')
        #Project each point and check result for correctness.
        CoordTransformer = qgis.core.QgsCoordinateTransform(sourceCRS, targetCRS, qgis.core.QgsProject.instance())
        self.assertTrue(CoordTransformer.isValid())
        self.assertFalse(CoordTransformer.isShortCircuited())
        sCRS_points=[]
        tCRS_points=[]
        features  = pointVectorLayer.getFeatures()
        for feature in features:
            sCRS_point = feature.geometry().asPoint()
            sCRS_points.append(sCRS_point)
            tCRS_points.append(CoordTransformer.transform(sCRS_point))
        #Check each projected point for correctness.
        for sCRS_point, tCRS_point in zip(sCRS_points, tCRS_points):
            longitude = sCRS_point.x()
            latitude = sCRS_point.y()
            expected_ksi = -np.tan(np.radians(-latitude/2 + 45))*np.cos(np.radians(longitude))
            expected_eta = -np.tan(np.radians(-latitude/2 + 45))*np.sin(np.radians(longitude))
            ksi = tCRS_point.x()
            eta = tCRS_point.y()
            ksiError = abs(expected_ksi - ksi)
            etaError = abs(expected_eta - eta)
            self.assertTrue(-1e-9 < ksiError < 1e-9)
            self.assertTrue(-1e-9 < etaError < 1e-9)


# these are a bunch of other WKT strings that might be useful in future.
        #targetCRS.createFromWkt('PROJCS["unnamed",GEOGCS["Normal Sphere (r=0.5)",DATUM["unknown",SPHEROID["sphere",6.37101e+6,0]],PRIMEM["Greenwich",0],UNIT["degree",0.01745329251994328]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",90],PARAMETER["central_meridian",90],UNIT["Meter",1]]')
        #Project each point and check result for correctness.
        #targetCRS.createFromWkt('PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6.37101e+6,0,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",90],PARAMETER["central_meridian",90],UNIT["Meter",1]]')
        #targetCRS.createFromWkt('PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",90],PARAMETER["central_meridian",90],UNIT["Meter",1]]')

        #sourceCRS_WKT = 'GEOGCS["Normal Sphere (r=6.37101e+6)",DATUM["WGS 84",SPHEROID["sphere",6.37101e+6,0]],PRIMEM["Greenwich",0],UNIT["degree",0.01745329251994328]]'
        #sourceCRS_WKT = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6.37101e+6,0,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

