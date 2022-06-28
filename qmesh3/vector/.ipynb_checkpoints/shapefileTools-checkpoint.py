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

from __future__ import absolute_import
import numpy as np
import os
import sys
from ..lib import Trees as Trees
import qgis.core


import logging
LOG = logging.getLogger(__package__)

class Shapes(object):
    ''' Wrapper object for ESRI shapefiles. '''
    def __init__(self,
                 coordRefSystem=None,
                 southPoleCoordinates=None,
                 isGlobal=False,
                 shapeType=None,
                 fields=qgis.core.QgsFields(),
                 featureList=[]):
        """ Initialise the Shapes object. """
        self.setCoordRefSystem(coordRefSystem)
        self.southPoleCoordinates = southPoleCoordinates
        self.isGlobal = isGlobal
        self.shapeType = shapeType
        self.fields = fields
        self.featureList = []
        self.featureCount = len(featureList)
        # Wrap around a vector layer.
        self.layer = None
        
    def setCoordRefSystem(self, coordRefSystem):
        """ Set the coordinate reference system from a QGIS CRS object.
        
        :arg qgis.core.QgsCoordinateReferenceSystem coordRefSystem: The coordinate reference system from QGIS.
        :returns: None 
        """
        # Check coordRefSystem is of type QgsCoordinateReferenceSystem
        # or None, then initialise variable.
        if coordRefSystem != None and  \
           type(coordRefSystem) != qgis.core.QgsCoordinateReferenceSystem:
            msg = 'Error: variable coordRefSystem is of incorrect type,'+\
                  str(type(coordRefSystem))
            raise BadArguments(msg)
        self.coordRefSystem = coordRefSystem
        
    def setCoordRefSystemFromString(self, coordRefSystem_string):
        """ Set the coordinate reference system based on the CRS' name. 
        
        :arg str coordRefSystem_string: The name of the CRS.
        :returns: None
        """
        if type(coordRefSystem_string) != str:
            raise BadArguments('Error: Input agrument coordRefSystem_string is of'+\
                               ' incorrect type, expected a string.\n')
        else:
            self.coordRefSystem = qgis.core.QgsCoordinateReferenceSystem()
            self.coordRefSystem.createFromString(coordRefSystem_string)
            if not self.coordRefSystem.isValid():
                msg = 'Error: could not initialise coordinate'+\
                      ' reference system '+coordRefSystem_string+' .\n'
                raise BadArguments(msg)

    def getCoordRefSystem(self):
        """ Get the coordinate reference system that has been set. 
        
        :rtype: qgis.core.QgsCoordinateReferenceSystem
        :returns: The coordinate reference system currently in use.
        """
        return self.coordRefSystem
    
    def setShapeType(self, shapeType):
        """ Sets the shape-type in Well-Known-Binary form. 
        :arg shapeType : shapeType The shape-type in well-known-binary form
        """
        #Check that shapeType describes an addmissible type.
        expectedShapeTypes = [qgis.core.QgsWkbTypes.LineString,
                              qgis.core.QgsWkbTypes.LineString25D,
                              qgis.core.QgsWkbTypes.MultiLineString,
                              qgis.core.QgsWkbTypes.MultiLineString25D,
                              qgis.core.QgsWkbTypes.MultiPoint,
                              qgis.core.QgsWkbTypes.MultiPoint25D,
                              qgis.core.QgsWkbTypes.MultiPolygon,
                              qgis.core.QgsWkbTypes.MultiPolygon25D,
                              qgis.core.QgsWkbTypes.Point,
                              qgis.core.QgsWkbTypes.Point25D,
                              qgis.core.QgsWkbTypes.Polygon,
                              qgis.core.QgsWkbTypes.Polygon25D]
        if shapeType not in expectedShapeTypes:
            msg = 'Error: Argument shapeType does not describe an'+\
                  ' expected shape-type.'
            raise BadArguments(msg)
        self.shapeType = shapeType
        
    def getShapeType(self):
        """ Return the shape-type in Well-Known-Binary form.
        :returns: The shape-type of the current object in well-known-binary form"""
        return self.shapeType
        
    def setFields(self, fields):
        ''' Sets the fields of the shapes.
        :arg qgis.core.QgsFields fields: fields-object containing field names and types.'''
        if len(self.featureList)>0:
            msg = 'Error: fields cannot be set once the Shapes'+\
                  ' object has been populated with features.'
            raise Exception(msg)
        if type(fields) != qgis.core.QgsFields:
            msg = 'Error: fields argument is of incorrect type.'
            raise BadArguments(msg)
        self.fields = fields
        
    def getFields(self):
        ''' Returns the fields of the shapes.
        :rtype: qgis.core.QgsFields
        :returns: A fields-object containing field names and types.'''
        return self.fields
        
    def setFeatureList(self, featureList):
        ''' Sets all features in object from given list.
        :arg list featureList: A list of features, each feature must be of type qgis.core.QgsFeature'''
        self.featureList = []
        for feature in featureList:
            self.addFeature(feature)
            
    def getFeatures(self):
        ''' Returns a list of all features.
        :rtype: list
        :returns: A list of features, each feature must be of type qgis.core.QgsFeature'''
        return self.featureList
        
    def addFeature(self, feature):
        ''' Add a single feature. The feature must have the same shape-type as well
        as the same fields with the shapes-object it is appended to. Also the
        shapes-object must have a valid Coordinate Reference System initialised.
        :arg qgis.core.QgsFeature feature: The feature to be added.'''
        #Check shape-type & CRS are set: Check the object is
        # initialised.
        if self.coordRefSystem == None:
            msg = 'Error: Cannot populate object if a Coordinate'+\
                  ' Reference System has not been set.'
            raise Exception(msg)
        #Check input variable is of correct type
        if type(feature) != qgis.core.QgsFeature:
            msg = 'Error: Argument feature is of incorrect type.'
            raise BadArguments(msg)
        #Check feature has no more or less fields than those stored in
        # this object's fields attribute.
        if len(feature.attributes()) != 0: #Required, feature.fields() returns NoneType otherwise
            for field in feature.fields():
                if field not in self.fields.toList():
                    msg = 'Error: Inserted feature does not have expected fields.'
                    raise Exception(msg)
        #Append feature to the list-store of the object.
        self.featureList.append(feature)
        self.featureCount += 1
        
    def removeFeatureWithAttribute(self, attributeName, attributeValue):
        '''Removes feature having with given attribute.
        :arg string attributeName : The name of the identifying attribute.
        :arg attributeValue: The value of the identifying attribute.'''
        featuresToRemove = []
        for feature in self.featureList:
            thisAttributeValue = feature.attribute(attributeName)
            if thisAttributeValue == attributeValue:
                featuresToRemove.append(feature)
        for feature in featuresToRemove:
            self.featureList.remove(feature)
                
    def removeFeature(self, feature):
        """ Remove a feature from the shapes object. 
        
        :arg feature: The feature to be removed. This must be present in the featuresList.        
        """
        # Check feature actually exists in object
        if feature not in self.featureList:
            msg = 'Attempting to remove a non-existing feature from a shapes object.'
            raise Exception(msg)
        self.featureList.remove(feature)
        
    def fromQgisLayer(self, shapefileVectorlayer):
        '''Set Coordinate Reference System, shape-type, fields and features from given
        qgis vector layer.
        :arg qgis.core.QgsVectorLayer shapefileVectorlayer: The layer to read from.'''
        import copy
        if not shapefileVectorlayer.isValid():
            msg = 'Error: Input vector layer is not valid.'
            raise Exception(msg)
        #Extract layer Coordinate Reference System. Then assign
        # object's property
        crs = shapefileVectorlayer.crs()
        self.setCoordRefSystem(crs)
        #Extract type of geometry in layer.
        self.setShapeType(shapefileVectorlayer.wkbType())
        #Extract the fields from the layer.
        self.layer = shapefileVectorlayer
        layerFields = shapefileVectorlayer.dataProvider().fields()
        self.setFields(layerFields)
        #Extract the features from the layer.
        features = shapefileVectorlayer.getFeatures()
        self.setFeatureList(features)
        
    def fromFile(self, inputFileName):
        '''Set Coordinate Reference System, shape-type, fields and features from given
        shapefile.
        :arg string inputFileName: Filename of shapefile to read from.'''
        #Check inputFileName is of type string, then
        # initialise variable.
        if type(inputFileName) != str:
            msg = 'Error: variable inputFileName is of incorrect type.'
            raise BadArguments(msg)
        #Check file exists
        if not os.path.isfile(inputFileName):
            msg = 'File '+inputFileName+' not found.'
            raise BadArguments(msg)
        #Open shapefile and create a vector-layer from it.
        LOG.info('Reading file '+inputFileName)
        shapefileVectorlayer = qgis.core.QgsVectorLayer(inputFileName,"Shapefile","ogr")
        self.fromQgisLayer(shapefileVectorlayer)
        
    def writeInfo(self):
        """ Display helpful debugging information about the CRS, the features, and geometry type. """
        
        # Display description of coordinate reference system.
        if self.coordRefSystem.description()!='' :
             LOG.debug('    Coordinate Reference System: '+\
                              self.coordRefSystem.description()+\
                              ' ('+str(self.coordRefSystem.authid())+').')
        elif self.coordRefSystem.toWkt()!='' :
             LOG.debug('Found coordinate reference'+
                              ' system with WKT: '+self.coordRefSystem.toWkt())
        else:
             LOG.debug('Found no coordinate reference system.')
        if self.coordRefSystem.description()!='':
            LOG.debug('Coordinate reference system has WKT: '+self.coordRefSystem.toWkt())
        # Print number of features.
        LOG.info('Number of features: '+str(self.featureCount))
        # Display type of geometry.
        if self.shapeType == qgis.core.QgsWkbTypes.LineString:
            LOG.debug(', of type line-string.')
        elif self.shapeType == qgis.core.QgsWkbTypes.LineString25D:
            LOG.debug(', of type 2.5 dim. line-string.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiLineString:
            LOG.debug(', of type multi-line-string.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiLineString25D:
            LOG.debug(', of type 2.5 dim. multi-line-string.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiPoint:
            LOG.debug(', of type multi-point.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiPoint25D:
            LOG.debug(', of type 2.5 dim. multi-point.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon:
            LOG.debug(', of type multi-polygon.')
        elif self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon25D:
            LOG.debug(', of type 2.5 dim. multi-polygon.')
        elif self.shapeType == qgis.core.QgsWkbTypes.Point:
            LOG.debug(', of type point.')
        elif self.shapeType == qgis.core.QgsWkbTypes.Point25D:
            LOG.debug(', of type 2.5 dim. point.')
        elif self.shapeType == qgis.core.QgsWkbTypes.Polygon:
            LOG.debug(', of type polygon.')
        elif self.shapeType == qgis.core.QgsWkbTypes.Polygon25D:
            LOG.debug(', of type 2.5 dim. polygon.')
        else:
            LOG.debug('.')
        # Print some information on each feature.
        for feature in self.featureList:
            # fetch geometry
            featGeometry = feature.geometry()
            LOG.debug('Feature '+str(feature.id())+ ' : ')
            # show some information about the feature
            if featGeometry.wkbType() == qgis.core.QgsWkbTypes.Point or \
               featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPoint or \
               featGeometry.wkbType() == qgis.core.QgsWkbTypes.Point25D or \
               featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPoint25D:
                if featGeometry.wkbType() ==  qgis.core.QgsWkbTypes.MultiPoint or \
                   featGeometry.wkbType() ==  qgis.core.QgsWkbTypes.MultiPoint25d:
                    point = featGeometry.asMultiPoint()
                else:
                    point = featGeometry.asPoint()
                LOG.debug('Point.')
                LOG.debug('Fields :')
                for field in feature.fields():
                    LOG.debug('            '+field.name()+' : ' + str(feature.attribute(field.name())))
                LOG.debug('                    Point Coords: '+point.toString())
            elif featGeometry.wkbType() == qgis.core.QgsWkbTypes.LineString or\
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString or \
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.LineString25D or \
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString25D:
                if featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString or \
                   featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString25D:
                    line = featGeometry.asMultiPolyline()
                else:
                    line = featGeometry.asPolyline()
                LOG.debug('Line with '+str(len(line))+' points.') 
                LOG.debug('            Fields :' )
                for field in feature.fields():
                   LOG.debug('            '+field.name()+' : ' + str(feature.attribute(field.name())))
            elif featGeometry.wkbType() == qgis.core.QgsWkbTypes.Polygon or \
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPolygon or \
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.Polygon25D or \
                 featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPolygon25D:
                if featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPolygon or \
                   featGeometry.wkbType() == qgis.core.QgsWkbTypes.MultiPolygon25D:
                    polygon = featGeometry.asMultiPolygon()
                else:
                    polygon = featGeometry.asPolygon()
                numPts = 0
                for ring in polygon:
                    numPts += len(ring)
                LOG.debug('        Polygon with '+str(len(polygon))+ ' rings and '+str(numPts)+ ' points.')
                LOG.debug('            Fields :' )
                for field in feature.fields():
                   LOG.debug('            '+field.name()+' : ' + str(feature.attribute(field.name())))
                for (ringIndex, ring) in zip(list(range(len(polygon))), polygon):
                   LOG.debug('                Ring '+str(ringIndex)+' :')
                   for point in ring:
                      LOG.debug('                Point Coords:'+point.toString())
            else:
                LOG.warning('        Unknown geometry type!')

    def writeFile(self, outputFileName):
        '''Write shapefile from current shapes. Overwrites existing files.
        :arg string outputFileName: Filename of shapefile to write to.'''
        LOG.info('Writing file '+outputFileName+' ...')
        # Write output shapefile: Open file and check for writer error.
        # Then write features and close file.
        fileWriter = qgis.core.QgsVectorFileWriter(outputFileName,
                   "CP1250", self.fields, self.shapeType, self.getCoordRefSystem(),
                   "ESRI Shapefile")
        if fileWriter.hasError() != qgis.core.QgsVectorFileWriter.NoError:
            raise Exception('Error when creating shapefile '+outputFileName+
                        ' : '+str(fileWriter.hasError()))
        for feature in self.featureList:
            fileWriter.addFeature(feature)
        del fileWriter
        
    def getFeatureCount(self):
        """ Return the number of features currently in the Shapes object.
        
        :rtype: int
        :returns: The number of features.
        """
        return self.featureCount

    def getBoundingRectangleExtents(self):
        firstFeature = self.featureList[0]
        pointList = firstFeature.geometry().asPolyline()
        firstPoint = pointList[0]
        xMin = firstPoint.x()
        xMax = xMin
        yMin = firstPoint.y()
        yMax = yMin
        for feature in self.featureList:
            pointList = feature.geometry().asPolyline()
            for point in pointList:
                x = point.x()
                y = point.y()
                if x > xMax: xMax=x
                if x < xMin: xMin=x
                if y > yMax: yMax=y
                if y < yMin: yMin=y
        return [xMin, xMax, yMin, yMax]
        
    def projectLonLatToUnitDisk(self, surfaceRadius=6.37101e+6, southPoleCoordinates=None):
        ''' Performs azimuthal projection to data, suitable for global geometries when
        working with planet-centred cartesian system. The geometry must be in EPSG:4326.
        While it is possible to convert non-EPSG:4326 systems to EPSG:4326
        and proceed, this is not done and left to the user to do it sepatatelly (via the changeCoordRefSystem method)
        :arg surfaceRadius float: The surface radious of the sphere approximating the Earth's surface. Defaults to 6.37101e+6
        :arg southPoleCoordinates float: The Coordinates of the south pole. This should only be used when meshing geometries
        where no land-mass is present at the South pole; in this case the southPoleCoordinates argument should specify the coordinates
        of a point on a land-mass.'''
        if self.getCoordRefSystem().authid() != 'EPSG:4326':
            msg = 'ERROR: Global files can only be processed if the'+\
                  ' Coordinate Reference System is EPSG:4326.'
            raise Exception(msg)
        #Todo: set CRS to user-specified.
        for feature in self.featureList:
            pointList = feature.geometry().asPolyline()
            for point in pointList:
                #If a rotation is required, convert coordinates to
                # cartesian (assuming a perfect sphere of unit radius), apply rotation
                # and convert back to lon-lat.
                if southPoleCoordinates != None:
                    longitude = point.x()
                    latitude = point.y()
                    x = np.cos(np.radians(longitude))*np.cos(np.radians(latitude))
                    y = np.sin(np.radians(longitude))*np.cos(np.radians(latitude))
                    z = np.sin(np.radians(latitude))
                    [X,Y,Z] = rotateCartesianBasis([x,y,z], southPoleCoordinates)
                    longitude = np.degrees(np.arctan2(Y, X))
                    latitude = np.degrees(np.arcsin(Z))
                else:
                    longitude = point.x()
                    latitude = point.y()
                r = -(1./180.)*latitude + 0.5
                disk_x = r*np.cos(np.radians(longitude))
                disk_y = r*np.sin(np.radians(longitude))
                point.setX(disk_x)
                point.setY(disk_y)
            feature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointList))
            
    def projectUnitDiskToLonLat(self):
        '''Converts geometry into EPSG:4326 from an azimuthal projection. Used when meshing global geometries.'''
        #Ensure the shape-type of the features is single-type, not mutli-type
        self.decomposeMultiFeatures()
        #Todo: set CRS from user-specified to EPSG:4326.
        for feature in self.featureList:
            featureGeometry = feature.geometry()
            shapeType = featureGeometry.wkbType()
            if shapeType == qgis.core.QgsWkbTypes.Polygon or \
               shapeType == qgis.core.QgsWkbTypes.Polygon25D:
                polygon = featureGeometry.asPolygon()
                for pointList in polygon:
                    for point in pointList:
                        disk_x = point.x()
                        disk_y = point.y()
                        r = np.sqrt(np.power(disk_x,2) + np.power(disk_y,2))
                        latitude = (r - 0.5)*-180.0
                        longitude = np.degrees(np.arctan2(disk_y,disk_x))
                        point.setX(longitude)
                        point.setY(latitude)
                feature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY(polygon))
            elif shapeType == qgis.core.QgsWkbTypes.Point or \
                 shapeType == qgis.core.QgsWkbTypes.Point25D:
                pointList = featureGeometry.asPoint()
                for point in pointList:
                    disk_x = point.x()
                    disk_y = point.y()
                    r = np.sqrt(np.power(disk_x,2) + np.power(disk_y,2))
                    latitude = (r - 0.5)*-180.0
                    longitude = np.degrees(np.arctan2(disk_y,disk_x))
                    point.setX(longitude)
                    point.setY(latitude)
                feature.setGeometry(qgis.core.QgsGeometry.fromPolyline(pointList))
            elif shapeType == qgis.core.QgsWkbTypes.LineString or \
                 shapeType == qgis.core.QgsWkbTypes.LineString25D:
                pointList = featureGeometry.asPolyline()
                for point in pointList:
                    disk_x = point.x()
                    disk_y = point.y()
                    r = np.sqrt(np.power(disk_x,2) + np.power(disk_y,2))
                    latitude = (r - 0.5)*-180.0
                    longitude = np.degrees(np.arctan2(disk_y,disk_x))
                    point.setX(longitude)
                    point.setY(latitude)
                feature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointList))
                
    def changeCoordRefSystem(self, targetCRS_string):
        """ Change to a different coordinate reference system.
        
        :arg str targetCRS_string: The name of the desired coordinate reference system.
        """
        # Make sure input is correct.
        assert isinstance(targetCRS_string, str)
        targetCRS = qgis.core.QgsCoordinateReferenceSystem()
        targetCRS.createFromString(targetCRS_string)
        CoordTransformer = qgis.core.QgsCoordinateTransform(
                self.getCoordRefSystem(), targetCRS, qgis.core.QgsProject.instance() )
        # Check if coordinate transformation is correctly initialised
        if not(CoordTransformer.isValid()):
            LOG.error('Change of Coordinate Reference System failed.'+
                    ' Initialisation of coordinate transformation failed.')
            raise Exception('Error in changeCoordRefSystem. QgsCoordinateTransform is not valid')
        #Figure out shape-type
        if self.shapeType == qgis.core.QgsWkbTypes.MultiPoint or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPoint25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Point25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Point:
            shapeTypeString = 'point'
        elif self.shapeType == qgis.core.QgsWkbTypes.LineString or\
           self.shapeType == qgis.core.QgsWkbTypes.LineString25D or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiLineString or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiLineString25D:
            shapeTypeString = 'line'
        elif self.shapeType == qgis.core.QgsWkbTypes.Polygon or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Polygon25D:
            shapeTypeString = 'polygon'
        # Ensure the shape-type of the features is single-type, not mutli-type
        self.decomposeMultiFeatures()
        # Transform point coordinates
        LOG.info('Converting '+shapeTypeString+' geometry to target reference coordinate system: '+targetCRS_string+'...')
        for feature in self.featureList:
            featureGeometry = feature.geometry()
            shapeType = featureGeometry.wkbType()
            # When the feature shape-type is a polygon, iterate through a list
            # of lists, where the first inner list encodes the outer boundary
            # and the following inner lists encode holes in the polygon.
            if shapeType == qgis.core.QgsWkbTypes.Polygon or \
               shapeType == qgis.core.QgsWkbTypes.Polygon25D:
                polygon = feature.geometry().asPolygon()
                for pointList in polygon:
                    for point in pointList:
                        newPoint = CoordTransformer.transform(point)
                        point.setX(newPoint.x())
                        point.setY(newPoint.y())
                feature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY(polygon))
            # When the feature shape-type is a point-set, iterate through a list of points.
            elif shapeType == qgis.core.QgsWkbTypes.Point or \
                 shapeType == qgis.core.QgsWkbTypes.Point25D:
                pointList = feature.geometry().asPoint()
                for point in pointList:
                    newPoint = CoordTransformer.transform(point)
                    point.setX(newPoint.x())
                    point.setY(newPoint.y())
                feature.setGeometry(qgis.core.QgsGeometry.fromPolyline(pointList))
            # When the feature shape-type is a line, iterate through a list of points.
            elif shapeType == qgis.core.QgsWkbTypes.LineString or \
                 shapeType == qgis.core.QgsWkbTypes.LineString25D:
                pointList = feature.geometry().asPolyline()
                for point in pointList:
                    newPoint = CoordTransformer.transform(point)
                    point.setX(newPoint.x())
                    point.setY(newPoint.y())
                feature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointList))
            else:
                attrs = feature.attributes()
                msg = 'Error in Coordinate Reference System Change: ' +\
                        'shape-type of object was not admitted: ' + str(shapeType) +\
                        ' ' + str(attrs)
                raise Exception(msg)
        # Assign CRS of transformed coordinates to object attribute.
        self.setCoordRefSystem(targetCRS)

    def decomposeMultiFeatures(self):
        '''Decompose multi-features into sinlge-features. Works on point, line and polygon type features'''
        if self.shapeType == qgis.core.QgsWkbTypes.MultiPoint or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPoint25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Point25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Point:
            shapeTypeString = 'point'
        elif self.shapeType == qgis.core.QgsWkbTypes.LineString or\
           self.shapeType == qgis.core.QgsWkbTypes.LineString25D or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiLineString or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiLineString25D:
            shapeTypeString = 'line'
        elif self.shapeType == qgis.core.QgsWkbTypes.Polygon or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon or\
           self.shapeType == qgis.core.QgsWkbTypes.MultiPolygon25D or\
           self.shapeType == qgis.core.QgsWkbTypes.Polygon25D:
            shapeTypeString = 'polygon'
        LOG.debug('        Breaking up '+shapeTypeString+' multi-features...')
        #Locate multi-features and decompose.
        fields = self.getFields()
        newFeatures = []
        multiFeatures = []
        singleFeaturesCount = 0
        for feature in self.featureList:
            featureGeometry = feature.geometry()
            shapeType = featureGeometry.wkbType()
            #Locate multi-points.
            if shapeType == qgis.core.QgsWkbTypes.MultiPoint or \
               shapeType == qgis.core.QgsWkbTypes.MultiPoint25D:
                multiPoints = featureGeometry.asMultiPoint()
                for point in multiPoints:
                    newFeature = qgis.core.QgsFeature()
                    newfeature.setGeometry(qgis.core.QgsGeometry.fromPoint(point))
                    #Set fields and attributes of new lines, by copying the fields
                    # and attributes of the multi-line.
                    newFeature.setFields(fields)
                    for field in feature.fields().toList():
                        attributeIdentifier = field.name()
                        attributeValue = feature.attribute(attributeIdentifier)
                        newFeature.setAttribute(attributeIdentifier, attributeValue)
                    newFeatures.append(newFeature)
                multiFeatures.append(feature)
            #Locate multi-lines.
            if shapeType == qgis.core.QgsWkbTypes.MultiLineString or \
               shapeType == qgis.core.QgsWkbTypes.MultiLineString25D:
                multiLines = featureGeometry.asMultiPolyline()
                for line in multiLines:
                    newFeature = qgis.core.QgsFeature()
                    newFeature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(line))
                    #Set fields and attributes of new lines, by copying the fields
                    # and attributes of the multi-line.
                    newFeature.setFields(fields)
                    for field in feature.fields().toList():
                        attributeIdentifier = field.name()
                        attributeValue = feature.attribute(attributeIdentifier)
                        newFeature.setAttribute(attributeIdentifier, attributeValue)
                    newFeatures.append(newFeature)
                multiFeatures.append(feature)
            #Locate multi-polygons.
            if shapeType == qgis.core.QgsWkbTypes.MultiPolygon or \
               shapeType == qgis.core.QgsWkbTypes.MultiPolygon25D:
                multiPolygons = featureGeometry.asMultiPolygon()
                for polygon in multiPolygons:
                    newFeature = qgis.core.QgsFeature()
                    newFeature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY(polygon))
                    #Set fields and attributes of new lines, by copying the fields
                    # and attributes of the multi-line.
                    newFeature.setFields(fields)
                    for field in feature.fields().toList():
                        attributeIdentifier = field.name()
                        attributeValue = feature.attribute(attributeIdentifier)
                        newFeature.setAttribute(attributeIdentifier, attributeValue)
                    newFeatures.append(newFeature)
                multiFeatures.append(feature)
            else:
                singleFeaturesCount += 1
        #Remove the multi-features from object.
        for feature in multiFeatures:
            self.removeFeature(feature)
        #Add the new lines to object store.
        for feature in newFeatures:
            self.addFeature(feature)
        #Inform user
        LOG.debug('            Found '+str(singleFeaturesCount)+' single-part features.')
        LOG.debug('            Found '+str(len(multiFeatures))+' multi-part features,')
        LOG.debug('                decomposed into '+str(len(newFeatures))+' single-part features.')
        LOG.debug('            Now '+str(self.getFeatureCount())+' single-part features in total.')

    def makeLinesContiguous(self):
        '''Makes lines connected: Locates lines whose ends are not connected to end-points of other,
        lines and move the disconnected end points to the closest end-point. When isGlobal flag is set
        the method will correctly connect lines across the whole globe.'''
        LOG.info('    Making lines contiguous...')
        #If file describes a whole globe, CRS must be EPSG:4326.
        # Then use an azimuthal projection to ensure we can connect
        # lines that cross the ends of the map.
        if self.isGlobal:
            LOG.debug('        Making contiguous lines: Input describes'+\
                             ' a whole globe...')
            if self.getCoordRefSystem().authid()!='EPSG:4326':
                msg = 'Global geometries must be given in EPSG:4326. The Coordinate '+\
                      ' Reference System detected is '+\
                      contiguousLinesShapes.getCoordRefSystem().authid()+' .'
                LOG.error(msg)
                raise Exception('Incorrect Coordinate Reference System.')
            self.projectLonLatToUnitDisk()
        #Separate lines into contiguous and non-contiguous. Here, a line
        # qualifies as contiguous if their end-points are end- points in other
        # lines, or the line is a loop. Otherwise they are put into the
        # non-contiguous list.
        contiguous_lineFeatures = []
        non_contiguous_lineFeatures = {}
        for currentLine in self.featureList:
            currentLineGeometry = currentLine.geometry()
            shapeType = currentLineGeometry.wkbType()
            #Make sure the method is invoked on single-part line shapes
            if not(shapeType != qgis.core.QgsWkbTypes.LineString or \
                   shapeType != qgis.core.QgsWkbTypes.LineString25D):
                msg = 'Method makeLinesContiguous can only be invoked on single-part line objects.'
                raise Exception(msg)
            #Get the end-points of the current line
            pointsList = currentLine.geometry().asPolyline()
            startPoint = pointsList[0]
            endPoint = pointsList[-1]
            foundStartPoint=False
            foundEndPoint=False
            #If the starting point of the current line is the same as the
            # end-point, then we have identified a loop composed of a
            # single line.
            if startPoint == endPoint:
                foundStartPoint=True
                foundEndPoint=True
            #Otherwise, look into the remaining lines to find
            # connected lines. Note that the start-point and
            # end-point of the current line could connect to 
            # different lines.
            else:
                for otherLine in self.featureList:
                    if otherLine != currentLine:
                        pointsList = otherLine.geometry().asPolyline()
                        otherLine_startPoint = pointsList[0]
                        otherLine_endPoint = pointsList[-1]
                        if startPoint in [otherLine_startPoint, otherLine_endPoint]:
                            foundStartPoint=True
                        if endPoint in [otherLine_startPoint, otherLine_endPoint]:
                            foundEndPoint=True
            if foundStartPoint and foundEndPoint:
                contiguous_lineFeatures.append(currentLine)
            else:
                #Non-contiguous lines are stored in a map, alongisde flags indicating
                # which of their end-points are not connected. This way when we later
                # try to connect disconnected lines, we do not have to re-connect the
                # already connected ends.
                non_contiguous_lineFeatures[currentLine] = [foundStartPoint, foundEndPoint]
        LOG.debug('        Making contiguous lines: '+\
                           str(len(self.featureList)) +\
                         ' lines in total.')
        LOG.debug('        Found '+\
                           str(len(contiguous_lineFeatures)) +\
                         ' lines that are parts of contiguous segments.')
        LOG.debug('        Found '+str(len(non_contiguous_lineFeatures)) +\
                         ' lines not connected to other lines.')

        #Iterate though non-contiguous lines and connect their end-points
        # to other lines. Use a proximity criterion to connect end-points.
        LOG.debug('        Making contiguous lines: Connecting disconnected lines...')
        for currentLine in non_contiguous_lineFeatures:
            [foundStartPoint, foundEndPoint] = non_contiguous_lineFeatures[currentLine]
            pointsList = currentLine.geometry().asPolyline()
            #If the start point is disconnected, locate the closest point, by
            # searching in all lines. The search includes the current line, in the case
            # the current line was indended to be a 'loop'.
            if not foundStartPoint:
                #Get the start-point of the current line
                startPoint = pointsList[0]
                #Compute the distances between start-point of the current line
                # and end-points of the other line, while accumulating minima
                closestToStartPoint = None
                closestToStartPoint_distance = np.nan
                for otherLine in self.featureList:
                    #Get end-points of 'other' line
                    otherLine_pointsList = otherLine.geometry().asPolyline()
                    other_startPoint = otherLine_pointsList[0]
                    other_endPoint = otherLine_pointsList[-1]
                    #Distance between current-line start point and other-line start point
                    start_to_start_distance = \
                          np.sqrt(np.sum(np.power([startPoint.x() - other_startPoint.x(),
                                                   startPoint.y() - other_startPoint.y()],2.0)))
                    #Distance between current-line start point and other-line end point
                    start_to_end_distance = \
                          np.sqrt(np.sum(np.power([startPoint.x() - other_endPoint.x(),
                                                   startPoint.y() - other_endPoint.y()],2.0)))
                    #Compare distances and store closest point, carefull with the case
                    # of other line being the same as the current line.
                    if otherLine == currentLine:
                        if closestToStartPoint_distance > start_to_end_distance or \
                           np.isnan(closestToStartPoint_distance):
                            closestToStartPoint = other_endPoint
                            closestToStartPoint_distance = start_to_end_distance
                    else:
                        if start_to_start_distance < start_to_end_distance and \
                           closestToStartPoint == None:
                            closestToStartPoint = other_startPoint
                            closestToStartPoint_distance = start_to_start_distance
                        elif start_to_start_distance > start_to_end_distance and \
                           closestToStartPoint == None:
                            closestToStartPoint = other_endPoint
                            closestToStartPoint_distance = start_to_end_distance
                        elif start_to_start_distance < start_to_end_distance and \
                           closestToStartPoint_distance > start_to_start_distance:
                            closestToStartPoint = other_startPoint
                            closestToStartPoint_distance = start_to_start_distance
                        elif start_to_start_distance > start_to_end_distance and \
                           closestToStartPoint_distance > start_to_end_distance:
                            closestToStartPoint = other_endPoint
                            closestToStartPoint_distance = start_to_end_distance
                #Connect the start-point of the current line by moving it onto the
                # closest point. Then reconstruct geometry definition
                startPoint.setX(closestToStartPoint.x())
                startPoint.setY(closestToStartPoint.y())
                currentLine.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointsList))
            if not foundEndPoint:
                #Get the end-point of the current line
                endPoint = pointsList[-1]
                #Compute the distances between end-point of the current line
                # and end-points of the other line, while accumulating minima
                closestToEndPoint = None
                closestToEndPoint_distance = np.nan
                for otherLine in self.featureList:
                    #Get end-points of 'other' line
                    otherLine_pointsList = otherLine.geometry().asPolyline()
                    other_startPoint = otherLine_pointsList[0]
                    other_endPoint = otherLine_pointsList[-1]
                    #Distance between current-line end point and other-line start point
                    end_to_start_distance = \
                          np.sqrt(np.sum(np.power([endPoint.x() - other_startPoint.x(),
                                                   endPoint.y() - other_startPoint.y()],2.0)))
                    #Distance between current-line end point and other-line end point
                    end_to_end_distance = \
                          np.sqrt(np.sum(np.power([endPoint.x() - other_endPoint.x(),
                                                   endPoint.y() - other_endPoint.y()],2.0)))
                    #Compare distances and store closest point, carefull with the case
                    # of other line being the same as the current line.
                    if otherLine == currentLine:
                        if closestToEndPoint_distance > end_to_start_distance or \
                           np.isnan(closestToEndPoint_distance):
                            closestToEndPoint = other_startPoint
                            closestToEndPoint_distance = end_to_start_distance
                    else:
                        if end_to_start_distance < end_to_end_distance and \
                           closestToEndPoint == None:
                            closestToEndPoint = other_startPoint
                            closestToEndPoint_distance = end_to_start_distance
                        elif end_to_start_distance > end_to_end_distance and \
                           closestToEndPoint == None:
                            closestToEndPoint = other_endPoint
                            closestToEndPoint_distance = end_to_end_distance
                        elif end_to_start_distance < end_to_end_distance and \
                           closestToEndPoint_distance > end_to_start_distance:
                            closestToEndPoint = other_startPoint
                            closestToEndPoint_distance = end_to_start_distance
                        elif end_to_start_distance > end_to_end_distance and \
                           closestToEndPoint_distance > end_to_end_distance:
                            closestToEndPoint = other_endPoint
                            closestToEndPoint_distance = end_to_end_distance
                #Connect the end-point of the current line by moving it onto the
                # closest point. Then reconstruct geometry definition
                endPoint.setX(closestToEndPoint.x())
                endPoint.setY(closestToEndPoint.y())
                currentLine.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointsList))
        #If global project back onto EPSG:4326
        if self.isGlobal:
            self.projectUnitDiskToLonLat()



def rotateCartesianBasis(positionVectorCartesian, southPoleCoordinates, invertRotation=False):
    '''Function rotating the planet-centered Cartesian basis of the input vector. The rotation angles are specified vie entering the coordinates of the point that will end on the south pole after the rotation.'''
    southPoleLon = southPoleCoordinates[0]
    southPoleLat = southPoleCoordinates[1]
    if not invertRotation:
        angleAboutY = 90.0 + southPoleLat
        angleAboutZ = -1.0*southPoleLon
        rotatedCoords = rotateAboutY(rotateAboutZ(positionVectorCartesian, angleAboutZ), angleAboutY)
    if invertRotation:
        angleAboutY = -(90 + southPoleLat)
        angleAboutZ = southPoleLon
        rotatedCoords = rotateAboutZ(rotateAboutY(positionVectorCartesian, angleAboutY), angleAboutZ)
    return rotatedCoords.tolist()

def rotateAboutZ(positionVectorCartesian, angleInDegrees):
    '''Euler rotation about z-axis.
    Todo: unit-test'''
    import numpy as np
    angleInRad = np.radians(angleInDegrees)
    rotationMatrix = np.array([[ np.cos(angleInRad),  -np.sin(angleInRad),  0.0],
                               [ np.sin(angleInRad),  np.cos(angleInRad),   0.0],
                               [0.0,                  0.0,                  1.0]
                              ])  
    rotatedVectorCartesian = np.dot(rotationMatrix, np.array(positionVectorCartesian))
    return rotatedVectorCartesian

def rotateAboutY(positionVectorCartesian, angleInDegrees):
    '''Euler rotation about y-axis.
    Todo: unit-test'''
    import numpy as np
    angleInRad = np.radians(angleInDegrees)
    rotationMatrix = np.array([[ np.cos(angleInRad),  0.0,                  np.sin(angleInRad) ],
                               [0.0,                  1.0,                  0.0                ],
                               [-np.sin(angleInRad),  0.0,                  np.cos(angleInRad) ]
                              ])
    rotatedVectorCartesian = np.dot(rotationMatrix, np.array(positionVectorCartesian))
    return rotatedVectorCartesian

def rotateAboutX(positionVectorCartesian, angleInDegrees):
    '''Euler rotation about x-axis.
    Todo: unit-test'''
    import numpy as np
    angleInRad = np.radians(angleInDegrees)
    rotationMatrix = np.array([[1.0,                  0.0,                  0.0],
                               [0.0,                  np.cos(angleInRad),   -np.sin(angleInRad)],
                               [0.0,                  np.sin(angleInRad),   np.cos(angleInRad)]
                              ])
    rotatedVectorCartesian = np.dot(rotationMatrix, np.array(positionVectorCartesian))
    return rotatedVectorCartesian

def readShapefile(filename,):
    """ Read in a shape file to create a new Shapes object. 
    
    :arg str filename: The name of the shape file to read in.
    :rtype: list
    :returns: A list containing the coordinate reference system and the shape file's vector layer.
    """
    shapes = Shapes()
    shapes.fromFile(filename)
    shapes.writeInfo()
    crs = shapes.getCoordRefSystem()
    shapefileVectorlayer = shapes.layer
    return [crs, shapefileVectorlayer]

def lines2polygon(outerRingFeatures, internalRings):
    '''Constructs a polygon from given lines.
    :arg qgis.core.QgsFeature outerRingFeatures: A list of line-features composing the outer boundary of the polygon.
    :arg qgis.core.QgsFeature outerRingFeatures: A list of rings (list of lists) composing holes in the polygon.
    :rtype: qgis.core.QgsFeature
    :returns: A polygon-feature'''
    LOG.debug('        Processing outer ring, composed of '+str(len(outerRingFeatures))+' features...')
    #Start with outer ring
    outerRingLines=[]
    for feature in outerRingFeatures:
        geometry = feature.geometry()
        # Ensure feature describes a line.
        if (not (geometry.wkbType() == qgis.core.QgsWkbTypes.LineString or
                 geometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString)):
            raise BadArguments('Input geometry does not appear to be of line-type.')
        if geometry.wkbType() == qgis.core.QgsWkbTypes.LineString:
            temp_lines = geometry
            temp_lines.convertToMultiType()
            temp_lines = temp_lines.asMultiPolyline()
            polyLine = []
            for t in temp_lines:
                polyLine.extend(t)
            outerRingLines.append(polyLine)
        else:
            outerRingLines.append(geometry.asMultiPolyline())
    LOG.debug('            Found %d lines in features of outer ring.' % len(outerRingLines))
    #Make list of outer lines consistent, and form outer ring of
    # the polygon-feature
    outerRing = []
    line = outerRingLines.pop()
    for point in line:
        outerRing.append(qgis.core.QgsPointXY(point))
    maxIter = len(outerRingLines)
    iterCount = 0
    while len(outerRingLines)>0:
        currentPoint = outerRing[-1]
        #Locate next line-segment on the boundary by checking the
        # coordinates of first and last points of all lines
        for line in outerRingLines:
            if (len(line) == 0):
                outerRingLines.remove(line)
                continue
            firstPoint = line[0]
            lastPoint = line[-1]
            if firstPoint == currentPoint:
                for point in line[1:]:
                    outerRing.append(qgis.core.QgsPointXY(point))
                outerRingLines.remove(line)
                break
            elif lastPoint == currentPoint:
                line.reverse()
                for point in line[1:]:
                    outerRing.append(qgis.core.QgsPointXY(point))
                outerRingLines.remove(line)
                break
        if (iterCount > maxIter):
            LOG.warning("Didn't finish processing all lines, but will soldier on regardless")
            break
        iterCount += 1
    #Check first and last points of the ring are identical:
    # Check the ring is closed!
    if outerRing[0] != outerRing[-1]:
        raise BadGeometry('Outer ring of polygon does not appear to be closed. First and last points must be identical.')

    #Start populating the polygon
    polygonList = [outerRing[0]]

    #Process the "holes" into rings
    LOG.debug('        Processing '+str(len(internalRings))+' internal rings...')
    for ringFeatures in internalRings:
        LOG.debug('            Processing internal ring, composed of '+str(len(ringFeatures))+' feature(s):')
        ringLines = []
        for feature in ringFeatures:
            LOG.debug('                ')
            for fieldIndex in np.arange(len(feature.fields())):
                fieldName = feature.fields()[fieldIndex].name()
                LOG.debug(fieldName+':'+str(feature.attribute(fieldName)))
                if fieldIndex == len(feature.fields())-1:
                    LOG.debug('\n')
                else:
                    LOG.debug(' , ')
            geometry = feature.geometry()
            #Ensure feature describes a line.
            if not(geometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString or
                   geometry.wkbType() == qgis.core.QgsWkbTypes.LineString):
                raise BadArguments('Input geometry does not appear to be of line-type.')
            if geometry.wkbType() == qgis.core.QgsWkbTypes.MultiLineString:
                ringLines.append(geometry.asMultiPolyline())
            if geometry.wkbType() == qgis.core.QgsWkbTypes.LineString:
                ringLines.append(geometry.asPolyline())

        #Make list of outer lines consistent, and a ring of
        # the polygon-feature
        ring = []
        line = ringLines.pop()
        for point in line:
            ring.append(point)
        while len(ringLines)>0:
            currentPoint = ring[-1]
            #Locate next line-segment on the boundary by checking the
            # coordinates of first and last points of all lines
            for line in ringLines:
                firstPoint = line[0]
                lastPoint = line[-1]
                if firstPoint == currentPoint:
                    ringLines.remove(line)
                    for point in line[1:]:
                        ring.append(point)
                if lastPoint == currentPoint:
                    ringLines.remove(line)
                    line.reverse()
                    for point in line[1:]:
                        ring.append(point)
        #Check first and last points of the ring are identical:
        # Check the ring is closed!
        if ring[0] != ring[-1]:
            raise BadGeometry('Ring of polygon does not appear to be closed. First and last points must be identical.')
        #Add ring to the polygon
        polygonList.append(ring[0])

    #Create a polygon-feature from the polygon-list and return.
    polygonFeature = qgis.core.QgsFeature()
    polygonFeature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY([polygonList]))
    return polygonFeature

def insertRegions(receivingLines, receivingPolygons,
                  insertedLines, insertedPolygons):
    '''Insert polygons inside given polygons. This method is intended for use
    in cases where the user wants to construct a hole (not meshed) into a given
    polygon and then construct a polygon, fitting that hole. Ideal for insterting
    regions parameterising tidal turbines into a meshed area. All polygon must
    be identified in terms of both polygon and line shapefiles. Both sets of input
    files must have PhysID defined and the same Coordinate Reference System.
    :arg Shapes receivingLines:
    :arg Shapes receivingPolygons:
    :arg Shapes insertedLines:
    :arg Shapes insertedPolygons:
    :rtype: tuple
    returns pair of Shapes objects, lines and polygons
    '''
    try:
        import PyQt4 as PyQt
    except ModuleNotFoundError:
        import PyQt5 as PyQt
    #Get Coordinate reference systems from the receiving and inserted 
    # shapes and check that the two sets of input shapes have the
    # same projection/CRS
    receiving_lines_crs = receivingLines.getCoordRefSystem()
    receiving_polygons_crs = receivingPolygons.getCoordRefSystem()
    inserted_lines_crs = insertedLines.getCoordRefSystem()
    inserted_polygons_crs = insertedPolygons.getCoordRefSystem()
    if receiving_lines_crs != receiving_polygons_crs or \
       receiving_lines_crs != inserted_lines_crs or \
       receiving_lines_crs != inserted_polygons_crs:
        raise Exception('Error: Input shapefiles do not have the same CRS')
    #Check features in all layers have at least the 'PhysID'
    # attribute defined
    for feature in receivingLines.getFeatures():
        #Create a list of the names of all the fields
        fieldNames = []
        for field in feature.fields():
            fieldNames.append(field.name())
        #Check 'PhysID' is in the above-created list.
        if 'PhysID' not in fieldNames:
            error_msg = 'Error: Feature in receiving lines '+\
                        ' does not have a PhysID attribute.\n'+\
                        '        Feature '+str(feature.id())+ '\n'
            for fieldName in fieldNames:
                error_msg += '            '+fieldName+' : ' +\
                       str(feature.attribute(fieldName))+'\n'
            raise Exception(error_msg)
    for feature in receivingPolygons.getFeatures():
        #Create a list of the names of all the fields
        fieldNames = []
        for field in feature.fields():
            fieldNames.append(field.name())
        #Check 'PhysID' is in the above-created list.
        if 'PhysID' not in fieldNames:
            error_msg = 'Error: Feature in receiving polygons '+\
                        ' does not have a PhysID attribute.\n'+\
                        '        Feature '+str(feature.id())+ '\n'
            for fieldName in fieldNames:
                error_msg += '            '+fieldName+' : ' +\
                       str(feature.attribute(fieldName))+'\n'
            raise Exception(error_msg)
    for feature in insertedLines.getFeatures():
        #Create a list of the names of all the fields
        fieldNames = []
        for field in feature.fields():
            fieldNames.append(field.name())
        #Check 'PhysID' is in the above-created list.
        if 'PhysID' not in fieldNames:
            error_msg = 'Error: Feature in inserted lines '+\
                        ' does not have a PhysID attribute.\n'+\
                        '        Feature '+str(feature.id())+ '\n'
            for fieldName in fieldNames:
                error_msg += '            '+fieldName+' : ' +\
                       str(feature.attribute(fieldName))+'\n'
            raise Exception(error_msg)
    for feature in insertedPolygons.getFeatures():
        #Create a list of the names of all the fields
        fieldNames = []
        for field in feature.fields():
            fieldNames.append(field.name())
        #Check 'PhysID' is in the above-created list.
        if 'PhysID' not in fieldNames:
            error_msg = 'Error: Feature in inserted polygons '+\
                        ' does not have a PhysID attribute.\n'+\
                        '        Feature '+str(feature.id())+ '\n'
            for fieldName in fieldNames:
                error_msg += '            '+fieldName+' : ' +\
                       str(feature.attribute(fieldName))+'\n'
            raise Exception(error_msg)
    #Line-insertion. Check inserted-lines do not intersect the
    # receiving-lines.
    for ins_feature in insertedLines.getFeatures():
        ins_featGeometry = ins_feature.geometry()
        for rec_feature in receivingLines.getFeatures():
            rec_featGeometry = rec_feature.geometry()
            if rec_featGeometry.intersects(ins_featGeometry):
                error_msg = 'Error: Receiving lines must not'+\
                            ' intersect inserted lines.\n'+\
                            'Intersecting features:\n'+\
                            '        Feature '+str(feature.id())+ '\n'
                for field in rec_feature.fields():
                    error_msg += '            '+field.name()+' : ' +\
                           str(rec_feature.attribute(field.name()))+'\n'
                error_msg += '        Feature '+str(feature.id())+ '\n'
                for field in ins_feature.fields():
                    error_msg += '            '+field.name()+' : ' +\
                           str(ins_feature.attribute(field.name()))+'\n'
                raise BadGeometry(error_msg)
    LOG.info('    Inserting '+str(insertedLines.getFeatureCount())+' line-features...')
    # Construct empty lists for output line-features and copy line features
    # from input objects, both the receiving and inserted lines.
    output_lineFeatures = []
    for rec_feature in receivingLines.getFeatures():
        output_lineFeatures.append(rec_feature)
    for ins_feature in insertedLines.getFeatures():
        output_lineFeatures.append(ins_feature)
    # Construct empty lists for output polygon-features. However, we
    # cannot just copy input polygon features, as receiving_polygonsLayer
    # and inserted_polygonsLayer may contain multiple polygons. We
    # must locate where each inserted polygon lies into, and insert it there.
    output_polygonFeatures = []
    for rec_feature in receivingPolygons.getFeatures():
        rec_featGeometry = rec_feature.geometry()
        for ins_feature in insertedPolygons.getFeatures():
            ins_featGeometry = ins_feature.geometry()
            rec_feature_counter = insertedPolygons.getFeatureCount()
            #If an inserted polygon is contained into one of
            # the receiving polygons, its outer ring is inserted as
            # a 'hole' to that receiving polygon. The inserted polygon
            # is also appended to the receiving feature list as a
            # separate polygon.
            if rec_featGeometry.intersects(ins_featGeometry):
                ins_featGeometry.convertToSingleType()
                ins_polygon = ins_featGeometry.asPolygon()
                ins_polygon_outerRing = ins_polygon[0]
                rec_featGeometry.addRing(ins_polygon_outerRing)
                # Re-assign geometry to receiving feature
                rec_feature.setGeometry(rec_featGeometry)
                #break
            else:
                rec_feature_counter -= 1
            #If an inserted polygon is not contained into one of
            # the receiving polygons, it is still inserted to the
            # output as a separate polygon.
            if rec_feature_counter == 0:
                stdout_msg = '                Following polygon does'+\
                             ' not appear to lie inside'+\
                             ' any of the receiving polygons.\n'+\
                             '                    Feature ' +\
                             str(feature.id())+ '\n'
                for field in ins_feature.fields():
                    stdout_msg += '                    '+\
                           field.name()+' : ' +\
                           str(ins_feature.attribute(field.name()))+'\n'
                stdout_msg += '               It will be inserted as a'+\
                              ' separate polygon.'
                LOG.debug(stdout_msg)
        # After making a "hole" to appropriate polygon, we can append polygon
        # into output polygon features list.
        output_polygonFeatures.append(rec_feature)
    # Append the polygons to fill the "holes" created above.
    for ins_feature in insertedPolygons.getFeatures():
        output_polygonFeatures.append(ins_feature)
    #Create Shapes-objects for output-lines and output-polygons.
    physicalIDfield = qgis.core.QgsField("PhysID",
                                    PyQt.QtCore.QVariant.Int)
    fields = receivingLines.getFields()
    outputLines = Shapes()
    outputLines.setCoordRefSystem(receiving_lines_crs)
    outputLines.setShapeType(qgis.core.QgsWkbTypes.LineString)
    outputLines.setFields(fields)
    outputLines.setFeatureList(output_lineFeatures)
    fields = receivingPolygons.getFields()
    outputPolygons = Shapes()
    outputPolygons.setCoordRefSystem(receiving_polygons_crs)
    outputPolygons.setShapeType(qgis.core.QgsWkbTypes.Polygon)
    outputPolygons.setFields(fields)
    outputPolygons.setFeatureList(output_polygonFeatures)
    return (outputLines, outputPolygons)

def identifyLoops(inputShapes,
                  isGlobal=False, defaultPhysID=None, fixOpenLoops=False,
                  extraPointsPerVertex=0):
    '''Method for automatic identification of loops, needed for automatic
    polygon detection (see method identifyPolygons) A loop is either the outer
    boundary of a polygon (areas to be meshed) or the boundaries of polygon holes.
    This method will first make sure all lines are connected and with then attempt
    sorting lines into loops. The method returns a Shapes object of line-type,
    containing a copy of the input line-features with the additional field "loopID",
    numbering the loop each line participates in. Note that each line is only allowed
    to participate in one loop. In cases where lines participate in multiple loops
    polygons must be defined manually, by the user.
    :arg Shapes inputShapes: Input shapes of line-type.
    :arg bool isGlobal: Flag for global domains, on-the-sphere
    :arg int defaultPhysID: Default physical ID, assigned to features without an existing physical ID 
    :arg bool fixOpenLoops: Flag for fixing disconnected lines.
    :arg int extraPointsPerVertex: Number of extra points ti insert, per line segment
    :rtype: Shapes
    :returns: A Shapes object containing a copy of input line features, with an extra attribute
    identifying the loop each line participates in.
    '''
    try:
        import PyQt4 as PyQt
    except ModuleNotFoundError:
        import PyQt5 as PyQt
    LOG.info('Locating line-loops...')
    #Make sure the method is invoked on line-type shapes
    if not(inputShapes.shapeType != qgis.core.QgsWkbTypes.LineString or \
           inputShapes.shapeType != qgis.core.QgsWkbTypes.LineString25D or \
           inputShapes.shapeType != qgis.core.QgsWkbTypes.MultiLineString or \
           inputShapes.shapeType != qgis.core.QgsWkbTypes.MultiLineString25D) :
        msg = 'Method identifyLoops can only be invoked on line-type objects.\n'
        raise Exception(msg)
    #Find out if "loopID" field exists in input shapes.
    if inputShapes.fields.indexFromName('LoopID') == -1:
        LOG.debug('    No loop-ID`s found in line-shapefile.')
        haveLoopID = False
    else:
        LoopIDIndex = inputShapes.fields.indexFromName('LoopID')
        LOG.debug('    Found loop-ID attribute in line-shapefile as field with index '+str(LoopIDIndex))
        haveLoopID = True
    #Find out if "physID" field exists
    if inputShapes.fields.indexFromName('PhysID') == -1:
        LOG.debug('    No physical ID`s found in line-shapefile.')
        haveLinePhysID = False
    else:
        PhysIDIndex = inputShapes.fields.indexFromName('PhysID')
        LOG.debug('    Found physical ID attribute in line-shapefile as field with index '+str(PhysIDIndex))
        haveLinePhysID = True
    #
    loopIDfield = qgis.core.QgsField("LoopID", PyQt.QtCore.QVariant.Int)
    physIDfield = qgis.core.QgsField("PhysID", PyQt.QtCore.QVariant.Int)
    fields = qgis.core.QgsFields()
    fields.append(loopIDfield)
    fields.append(physIDfield)
    #Exract line-features and store them in a list. First, note that
    # we are after a list of line features - not multi-line features.
    # Second, we need to set the fields to include both loopID and
    # physID, while the input might include just one or both. Thus below
    # we iterate and create new features, from the geometry of the
    # input.
    lineFeaturesList = []
    for feature in inputShapes.featureList:
        geom = feature.geometry()
        #Ensure feature describes a line. If a polyline, break down
        # into lines.
        if not(geom.wkbType() != qgis.core.QgsWkbTypes.LineString or \
               geom.wkbType() != qgis.core.QgsWkbTypes.MultiLineString):
            msg = 'Method identifyLoops can only be invoked on line-type'+\
                  ' objects.\n'
            raise Exception(msg)
        if geom.wkbType() == qgis.core.QgsWkbTypes.MultiLineString:
            temp_lines = geom.asMultiPolyline()
            for t in temp_lines:
                t = geom.fromPolylineXY(t)
                pointList = t.asPolyline()
                #Check if extra points per vertex are required, and insert
                if extraPointsPerVertex != 0:
                    pointListDense = []
                    numbOriginalPoints=len(pointList)
                    for pointCounter in range(numbOriginalPoints-1):
                        endPoint1 = pointList[pointCounter]
                        endPoint2 = pointList[pointCounter+1]
                        deltaXi = endPoint2.x() - endPoint1.x()
                        deltaEta = endPoint2.y() - endPoint1.y()
                        newNumbPointsPerVertex = extraPointsPerVertex + 2
                        for newPointCounter in range(newNumbPointsPerVertex):
                            newXi = endPoint1.x() + \
                                    (float(newPointCounter)/(float(newNumbPointsPerVertex)-1.))*deltaXi
                            newEta = endPoint1.y() + \
                                    (float(newPointCounter)/(float(newNumbPointsPerVertex)-1.))*deltaEta
                            newPoint = qgis.core.QgsPoint(newXi, newEta)
                            pointListDense.append(newPoint)
                    pointList = pointListDense
                newFeature = qgis.core.QgsFeature()
                try:
                    newFeature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointList))
                except TypeError:
                    newFeature.setGeometry(qgis.core.QgsGeometry.fromPolyline(pointList))
                newFeature.setFields(fields)
                if haveLoopID:
                    loopID = feature.attribute('LoopID')
                    newFeature.setAttribute('LoopID',loopID)
                if haveLinePhysID:
                    physID = feature.attribute('PhysID')
                    newFeature.setAttribute('PhysID',physID)
                lineFeaturesList.append(newFeature)
        else:
            pointList = feature.geometry().asPolyline()
            #Check if extra points per vertex are required, and insert
            if extraPointsPerVertex != 0:
                pointListDense = []
                numbOriginalPoints=len(pointList)
                for pointCounter in range(numbOriginalPoints-1):
                    endPoint1 = pointList[pointCounter]
                    endPoint2 = pointList[pointCounter+1]
                    deltaXi = endPoint2.x() - endPoint1.x()
                    deltaEta = endPoint2.y() - endPoint1.y()
                    newNumbPointsPerVertex = extraPointsPerVertex + 2
                    for newPointCounter in range(newNumbPointsPerVertex):
                        newXi = endPoint1.x() + \
                                (float(newPointCounter)/float(newNumbPointsPerVertex-1))*deltaXi
                        newEta = endPoint1.y() + \
                                (float(newPointCounter)/float(newNumbPointsPerVertex-1))*deltaEta
                        newPoint = qgis.core.QgsPoint(newXi, newEta)
                        pointListDense.append(newPoint)
                pointList = pointListDense
            newFeature = qgis.core.QgsFeature()
            newFeature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointList))
            newFeature.setFields(fields)
            if haveLoopID:
                loopID = feature.attribute('LoopID')
                newFeature.setAttribute('LoopID',loopID)
            if haveLinePhysID:
                physID = feature.attribute('PhysID')
                newFeature.setAttribute('PhysID',physID)
            lineFeaturesList.append(newFeature)
    #Assign PhysID attribute to each line. If attribute already exists
    # and a default-value is given, check for a null value and assign
    # the default-value instead.
    for feature in lineFeaturesList:
        if haveLinePhysID and defaultPhysID != None:
            if str(feature.attribute('PhysID')) == 'NULL':
                feature.setAttribute('PhysID', defaultPhysID)
        if not haveLinePhysID and defaultPhysID != None:
            feature.setAttribute('PhysID', defaultPhysID)
    #Populate a shapes object
    lineShapes = Shapes()
    lineShapes.setCoordRefSystem(inputShapes.coordRefSystem)
    lineShapes.setShapeType(qgis.core.QgsWkbTypes.LineString)
    lineShapes.setFields(fields)
    lineShapes.setFeatureList(lineFeaturesList)
    if isGlobal:
        lineShapes.isGlobal = isGlobal
    #If user wants to force-close open loops make the lines
    # contiguous
    if fixOpenLoops:
        lineShapes.makeLinesContiguous()
    #Locate line loops.
    #If file describes a whole globe, loop detection is easier
    # on a polar azimuthal equidistant projection. Also, input
    # CRS must be EPSG:4326. 
    if isGlobal:
        lineShapes.projectLonLatToUnitDisk()
    loopID = 0
    output_lineFeatures = []
    nonLoops = []
    while len(lineShapes.featureList) > 0:
        startingLine = lineShapes.featureList.pop()
        contiguousLines = []
        isLoop = False
        startingLine.setAttribute('LoopID',loopID)
        contiguousLines.append(startingLine)
        #Get an end-point of the current line
        pointsList = startingLine.geometry().asPolyline()
        startingLine_startPoint = pointsList[0]
        currentEndPoint = pointsList[-1]
        #If the starting point of the current line is the same as the
        # end-point, then we have identified a loop composed of a
        # single line.
        if startingLine_startPoint == currentEndPoint:
            isLoop = True
            output_lineFeatures.extend(contiguousLines)
            loopID += 1
        #Otherwise, look into the remaining lines to find out
        # connecting lines. We also introduce a safety counter,
        # to ensure that the loops are assembled or an error
        # message is given without the program hanging inside an
        # infinite loop. This is necessary as currently duplicate
        # and overlying lines are not treated.
        else:
            while not isLoop:
                foundNextLine = False
                for nextLine in lineShapes.featureList:
                    pointsList = nextLine.geometry().asPolyline()
                    if pointsList[-1] == currentEndPoint:
                        foundNextLine = True
                        nextLine.setAttribute('LoopID',loopID)
                        lineShapes.featureList.remove(nextLine)
                        #Reverse the order of the points, to ensure a consistent
                        # direction heading round the loop. The line is then appended
                        # to a list of contiguous lines. This way,
                        # contiguousLines[0].geometry().asPolyline()[0] can be thought
                        # as the "first" point on the contiguous segment and 
                        # contiguousLines[-1].geometry().asPolyline()[-1] can be seen
                        # as the last. Obviously, in a loop the two points should have
                        # the same coordinates.
                        pointsList.reverse()
                        nextLine.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(pointsList))
                        #Append this line to list of contiguous lines so far.
                        contiguousLines.append(nextLine)
                        currentEndPoint = pointsList[-1]
                        #Once the next line in the loop has been
                        # identified, check if we have closed the loop.
                        if currentEndPoint == startingLine_startPoint:
                            isLoop = True
                            output_lineFeatures.extend(contiguousLines)
                            loopID += 1
                            break
                    elif pointsList[0] == currentEndPoint:
                        foundNextLine = True
                        nextLine.setAttribute('LoopID',loopID)
                        lineShapes.featureList.remove(nextLine)
                        #Append this line to list of contiguous lines so far.
                        contiguousLines.append(nextLine)
                        currentEndPoint = pointsList[-1]
                        #Once the next line in the loop has been
                        # identified, check if we have closed the loop.
                        if currentEndPoint == startingLine_startPoint:
                            isLoop = True
                            output_lineFeatures.extend(contiguousLines)
                            loopID += 1
                            break
                #If the next line has not been found by the loop above the geometry is
                # still problematic!
                if not foundNextLine:
                    msg = 'Could not assemble loop. Line-sorting is hanging'+\
                          ' at loop '+str(loopID)+', at point'+str(currentEndPoint)+\
                          '. Duplicate or overlying'+\
                          ' lines could be present in the geometry.'
                    raise Exception(msg)
    #Check no lines intersect
    if fixOpenLoops:
        LOG.debug('        Checking for intersecting lines...')
    checkedLines=[]
    for lineCounter1 in range(len(output_lineFeatures)):
        for lineCounter2 in range(len(output_lineFeatures)):
            if lineCounter1 == lineCounter2 or\
               lineCounter2 in checkedLines: 
                break
            line1 = output_lineFeatures[lineCounter1]
            line2 = output_lineFeatures[lineCounter1]
            line1LoopID = line1.attribute('LoopID')
            line2LoopID = line2.attribute('LoopID')
            if line1.geometry().intersects(line2.geometry()):
                msg = 'WARNING: Output shapes contain intersecting lines:\n'+\
                      '         Line '+str(lineCounter1)+' with LoopID '+\
                      str(line1LoopID)+\
                      '         intersects line '+ str(lineCounter2)+\
                      ' with LoopID '+str(line2LoopID)
            checkedLines.append(lineCounter1)
    #Done looping, just output now.
    LOG.debug('        Found '+str(loopID)+' line-loops in total.')
    #Construct output Shapes object.
    outputShapes = Shapes()
    outputShapes.setCoordRefSystem(inputShapes.coordRefSystem)
    outputShapes.setShapeType(qgis.core.QgsWkbTypes.LineString)
    outputShapes.setFields(fields)
    outputShapes.setFeatureList(output_lineFeatures)
    #If shapes describe a whole globe, convert from unit-disk
    # to lon-lat EPSG:4326
    if isGlobal:
        outputShapes.projectUnitDiskToLonLat()
        outputShapes.isGlobal = isGlobal
    return outputShapes

def identifyPolygons(loopShapes,
                     isGlobal=False,
                     southPoleCoordinates=None,
                     meshedAreaPhysID = 0,
                     smallestNotMeshedArea=None,
                     smallestMeshedArea=None):
    '''Method for automatic polygon detection.
    :arg Shapes loopShapes: Input loop-shapes.
    :arg bool isGlobal: Flag for global domains, on-the-sphere.
    :arg tuple southPoleCoordinates: Location of point on a land-mass.
    :arg int meshedAreaPhysID: Physical ID attached to output polygons.
    :arg float smallestNotMeshedArea: Size of smallest hole.
    :arg float smallestMeshedArea: Size of smallest polygon.
    '''
    try:
        import PyQt4 as PyQt
    except ModuleNotFoundError:
        import PyQt5 as PyQt
    LOG.info('Building polygons...')
    #Ensure input Shapes object describes line-features.
    if loopShapes.shapeType != qgis.core.QgsWkbTypes.LineString:
        msg = 'Error, input shape-type must be QGis.WKBLineString.'
        raise BadGeometry(msg)
    #Ensure input Shapes object has a 'LoopID' field.
    if loopShapes.fields.indexFromName('LoopID') == -1:
        msg = 'Error, input Shapes-object does not have a "LoopID"'+\
              ' field.\n'
        raise BadGeometry(msg)
    #If shapes describe a whole globe, CRS must be EPSG:4326.
    # Check where the south-pole is and rotate coordinate system
    # if neccessary. Then, convert lon-lat to x-y on a unit-disk.
    if loopShapes.isGlobal:
        loopShapes.projectLonLatToUnitDisk(southPoleCoordinates=southPoleCoordinates)
        loopShapes.southPoleCoordinates = southPoleCoordinates
    #Extract features-list from input shapes, and buld a dictionary, as
    # loopID:[features of that loop]
    loops = {}
    inputFeatures = loopShapes.featureList
    for feature in inputFeatures:
        loopID = feature.attribute('LoopID')
        if loopID in list(loops.keys()):
            loops[loopID].append(feature)
        else:
            loops[loopID] = [feature]
    #Check for closed-ness of loops.
    for loopID in loops:
        loop = loops[loopID]
        startPoint = loop[0].geometry().asPolyline()[0]
        endPoint = loop[-1].geometry().asPolyline()[-1]
        if startPoint != endPoint:
            msg = 'ERROR: Input features with LoopID='+str(loopID)+\
                  ' does not form a closed loop, start & end points must'+\
                  ' be identical.\n'
            raise BadGeometry(msg)
    #For every loop, create a polygon.
    polygonDictionary={}
    for loopID in loops:
        loop = loops[loopID]
        #Create a single list containing all loop points. Recall a loop
        # is made up of multiple line-features
        loopPoints = []
        for line in loop:
            loopPoints.extend(line.geometry().asPolyline())
        #Check loop is composed of more than 2 points
        if len(loopPoints) <= 2:
            msg = 'Line-loop with LoopID '+str(loopID)+\
                  ' is composed of '+str(len(loopPoints))+' points.\n'+\
                  '         Will remove this line-loop and not construct '+\
                  'a polygon for this feature.'
            LOG.warning(msg)
            loopShapes.removeFeatureWithAttribute('LoopID', loopID)
            continue
        #Check that geometry resulting from points list is sane
        # do not add loop otherwise
        if qgis.core.QgsGeometry.fromPolygonXY([loopPoints]) == None:
            msg = 'Not possible to reconsitute a valid geometry'+\
                  ' from Line-loop with LoopID '+str(loopID)+'. \n'+\
                  '         Will remove this line-loop and not construct '+\
                  'a polygon for this feature.'
            LOG.warning(msg)
            loopShapes.removeFeatureWithAttribute('LoopID', loopID)
            continue
        #Passed all checks, add loop as polygon
        polygon = qgis.core.QgsFeature()
        polygon.setGeometry(qgis.core.QgsGeometry.fromPolygonXY([loopPoints]))
        polygonDictionary[loopID] = polygon
    #Locate outer-most polygons: They should not be contained by any
    # other polygons. Once located, store the polygons in new
    # dictionary, and remove them from polygonDictionary. This way,
    # we can iterate over polygonDictionary, locating the next set
    # of outer-most polygons that should describe islands. The third
    # iteration will recover lakes (water masses again). So on, until
    # the polygonDictionary becomes empty. In the meanwhile we built
    # a hierarchy of sea-masses - land-masses - sea-masses etc.
    polygonHierarchy = []
    while len(polygonDictionary)>0:
        #Locate all outer-most polygons.
        outerMostPolygons = {}
        for loopID in polygonDictionary:
            currentPolygon = polygonDictionary[loopID]
            isOuterPolygon = True
            for otherLoopID in polygonDictionary:
                if otherLoopID == loopID:
                    pass
                else:
                    otherPolygon = polygonDictionary[otherLoopID]
                    if otherPolygon.geometry().contains(currentPolygon.geometry()):
                        isOuterPolygon = False
                        break
            if isOuterPolygon:
                outerMostPolygons[loopID] = (currentPolygon)
        polygonHierarchy.append(outerMostPolygons)
        #Remove outer polygons from polygonDictionary
        for outerLoopID in outerMostPolygons:
            polygonDictionary.pop(outerLoopID)
    #Build polygon trees from the polygon hierarchy
    polygonTrees = []
    for hierarchyLevel, hierarchyLevelCounter in zip(polygonHierarchy, np.arange(len(polygonHierarchy))):
        #If no trees exist in list polygonTrees, Then we have just isolated the out-most boundaries of the domain.
        if len(polygonTrees) == 0:
            isMeshed = True
            for loopID, polygon in list(hierarchyLevel.items()):
                node = [loopID, polygon, isMeshed]
                polygonTrees.append(Trees.Tree(Trees.Branch((node))))
        #Otherwise locate the correct tree and add the polygon as a node.
        else:
            isParentMeshed = isMeshed
            isMeshed = not isParentMeshed
            for parentloopID, parentPolygon in list(polygonHierarchy[hierarchyLevelCounter-1].items()):
                for loopID, polygon in list(hierarchyLevel.items()):
                    if parentPolygon.geometry().contains(polygon.geometry()):
                        parentNode = [parentloopID, parentPolygon, isParentMeshed]
                        node = [loopID, polygon, isMeshed]
                        for tree in polygonTrees:
                            if tree.hasNode(parentNode):
                                tree.addNode(node, parentNode)
                                break

    #Inform user how many levels of nesting were found in the geometry
    LOG.debug('Number of outer boundaries found: '+str(len(polygonTrees)))
    #Initialise object used in subsequent area calculations
    areaCalculator = qgis.core.QgsDistanceArea()
    areaCalculator.setSourceCrs(loopShapes.getCoordRefSystem(),qgis.core.QgsCoordinateTransformContext())
    #areaCalculator.setSourceAuthId(loopShapes.getCoordRefSystem())
    areaCalculator.setEllipsoid(loopShapes.getCoordRefSystem().ellipsoidAcronym())
    #areaCalculator.setEllipsoidalMode(True)
    #areaCalculator.computeAreaInit()
    #Now the polygon trees are in place, insert the appropriate rings
    # to the correct polygons, to correctly define polygons to be
    # meshed/not-meshed. Also locate and remove islands that
    # fall below the user-specified threshold. Note when an island
    # is removed, all other polygons inside it are also removed from
    # the tree.
    removedNonMeshedCounter = 0
    #newPolygonTrees = []
    for tree in polygonTrees:
        for (node, childNodes) in tree.nodesWithChildren():
            parentPolygon = node[1]
            isParentMeshed = node[2]
            for childNode in childNodes:
                loopID = childNode[0]
                childPolygon = childNode[1]
                isChildMeshed = childNode[2]
                if smallestNotMeshedArea != None and not isChildMeshed:
                    area = areaCalculator.measureArea(childPolygon.geometry())
                    if area < smallestNotMeshedArea:
                        removed_subTree = tree.prune(childNode)
                        removedNonMeshedCounter += 1
                        #Corresponding entities must be removed from the loops too.
                        loopShapes.removeFeatureWithAttribute('loopID', loopID)
                    else:
                        new_feature = qgis.core.QgsFeature()
                        new_feature.setGeometry(parentPolygon.geometry())
                        new_feature_geom = new_feature.geometry()
                        new_feature_clipped = new_feature_geom.difference(childPolygon.geometry())
                        pp = qgis.core.QgsFeature()
                        pp.setGeometry(new_feature_clipped)
                        parentPolygon = pp
                        node[1] = parentPolygon
                else:
                    new_feature = qgis.core.QgsFeature()
                    new_feature.setGeometry(parentPolygon.geometry())
                    new_feature_geom = new_feature.geometry()
                    new_feature_clipped = new_feature_geom.difference(childPolygon.geometry())
                    pp = qgis.core.QgsFeature()
                    pp.setGeometry(new_feature_clipped)
                    parentPolygon = pp
                    node[1] = parentPolygon
            
    #Locate and remove water-masses that fall below the
    # user-specified threshold
    removedMeshedCounter = 0
    for tree in polygonTrees:
        for node in tree.nodes():
            polygon = node[1]
            isMeshed = node[2]
            if smallestMeshedArea != None and isMeshed:
                area = areaCalculator.measureArea(polygon.geometry())
                if area < smallestMeshedArea:
                    removed_subTree = tree.prune(node)
                    removedMeshedCounter += 1
                    #Corresponding entities must be removed from the loops too.
                    # Note that when removing lakes, we also need to remove
                    # any islands in that lake too. Tree-pruning above
                    # did this for the polygon tree. We must iterate through
                    # the removed sub-tree and remove all corresponding line-loops
                    # too.
                    for removedNode in removed_subTree.nodes():
                        loopID = removedNode[0]
                        isRemovedMeshed = removedNode[2]
                        loopShapes.removeFeatureWithAttribute('loopID', loopID)
                        if isRemovedMeshed:
                            removedMeshedCounter+=1
                        else:
                            removedNonMeshedCounter+=1
    # Removed lakes can create empty polygon trees, must now purge them.
    emptyTrees = []
    for tree in polygonTrees:
        if tree.isEmpty():
            emptyTrees.append(tree)
    for tree in emptyTrees:
            polygonTrees.remove(tree)
    LOG.info('Number of meshed-areas removed: '+str(removedMeshedCounter))
    LOG.info('Number of non-meshed areas removed: '+str(removedNonMeshedCounter))
    # Depending on verbocity, find meshed or non-meshed regions of
    # smallest area and report that metric back.
    # A log level of 10 corresponds to the DEBUG level.
    if LOG.getEffectiveLevel() == 10:
        meshedCounter = 0
        nonMeshedCounter = 0
        largestMeshedArea = None
        largestNotMeshedArea = None
        smallestMeshedArea = None
        smallestNotMeshedArea = None
        for tree in polygonTrees:
            for node in tree.nodes():
                polygon = node[1]
                isMeshed = node[2]
                area = areaCalculator.measureArea(polygon.geometry())
                #(area, units) = areaCalculator.convertMeasurement(
                #                       area,
                #                       qgis.core.QGis.Meters,
                #                       qgis.core.QGis.Degrees, True)
                #Leave area calculation to 
                # return in meters, for time being.
                if isMeshed:
                    meshedCounter += 1
                    if largestMeshedArea == None or largestMeshedArea < area:
                        largestMeshedArea = area
                    if smallestMeshedArea == None or smallestMeshedArea > area:
                        smallestMeshedArea = area
                else:
                    nonMeshedCounter += 1
                    if largestNotMeshedArea == None or largestNotMeshedArea < area:
                        largestNotMeshedArea = area
                    if smallestNotMeshedArea == None or smallestNotMeshedArea > area:
                        smallestNotMeshedArea = area
        LOG.debug('Number of areas to be meshed: '+str(meshedCounter))
        LOG.debug('Number of areas not to be meshed: '+str(nonMeshedCounter))
        LOG.debug('Largest meshed area: '+str(largestMeshedArea)+' m**2')
        LOG.debug('Smallest meshed area: '+str(smallestMeshedArea)+' m**2')
        LOG.debug('Largest not-meshed area: '+str(largestNotMeshedArea)+' m**2')
        LOG.debug('Smallest not-meshed area: '+str(smallestNotMeshedArea)+' m**2')
        
    # Populate a Shapes object with polygons to be meshed. Then return.
    outputShapes = Shapes()
    outputShapes.setCoordRefSystem(loopShapes.getCoordRefSystem())
    outputShapes.isGlobal = loopShapes.isGlobal
    outputShapes.southPoleCoordinates = loopShapes.southPoleCoordinates
    outputShapes.setShapeType(qgis.core.QgsWkbTypes.Polygon)
    physIDField = qgis.core.QgsField("PhysID",
                                    PyQt.QtCore.QVariant.Int)
    fields = qgis.core.QgsFields()
    fields.append(physIDField)
    outputShapes.setFields(fields)
    for tree in polygonTrees:
        for node in tree.nodes():
            polygon = node[1]
            isMeshed = node[2]
            polygon.setFields(fields)
            if isMeshed:
                polygon.setAttribute('PhysID',meshedAreaPhysID)
                outputShapes.addFeature(polygon)
                
    # If shapes describe a whole globe, convert from unit-disk to lon-lat EPSG:4326
    if isGlobal:
        outputShapes.projectUnitDiskToLonLat()
        loopShapes.projectUnitDiskToLonLat()
    return outputShapes

class BadArguments(Exception):
    def __init__(self, message):
        self.message=message
    def __str__(self):
        return self.message

class BadGeometry(Exception):
    def __init__(self, message):
        self.message=message
    def __str__(self):
        return self.message
