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

import numpy as np
import qgis.core

class Line(object):
    """ A class to represent line objects in the mesh. """
    def __init__(self, coordsList, coordRefSystem):
        """ Class initialisation method. 
        
        :arg list coordsList: The list of coordinates which make up the points in the line.
        :arg coordRefSystem: The coordinate reference system to use. This can be a qgis.core.QgsCoordinateReferenceSystem object, or a string representing the CRS' name.
        """
        if type(coordRefSystem) != str and \
           type(coordRefSystem) != qgis.core.QgsCoordinateReferenceSystem:
            raise BadArguments('Error: Input agrument coordRefSystem is of'+\
                               ' incorrect type, expected a string or'+\
                               ' QgsCoordinateReferenceSystem.\n')
        elif type(coordRefSystem) == str:
            self.coordRefSystem_string = coordRefSystem
            self.coordRefSystem = qgis.core.QgsCoordinateReferenceSystem()
            self.coordRefSystem.createFromString(coordRefSystem)
            if not self.coordRefSystem.isValid():
                msg = 'Error: could not initialise coordinate'+\
                      ' reference system '+self.coordRefSystem_string+' .\n'
                raise BadArguments(msg)
        else:
            self.coordRefSystem = coordRefSystem
            self.coordRefSystem_string = coordRefSystem.authid()
        #TODO: give a helpful message for empty coordsList
        self.coordsList = coordsList
        # Construct a QgsPoint-list from the coordinate-pairs. While the
        # same information is present in self.coordsList, we also
        # generate the pointsList, as it naturally goes together
        # with the coordinate reference system variable; especially
        # when changing the coordinate reference system. 
        self.pointsList = []
        for coordPair in self.coordsList:
            point = qgis.core.QgsPointXY()
            point.setX(coordPair[0])
            point.setY(coordPair[1])
            self.pointsList.append(point)
            
    def changeCoordRefSystem(self, targetCRS_string):
        """ Change the coordinate reference system.
        
        :arg str targetCRS_string: The name of the coordinate reference system to use.
        """
        raise Exception('Method implementation is incomplete.\n')
        
    def asCoordsList(self):
        """ Return the list of coordinates.
        
        :rtype: list
        :returns: The coordinates list.
        """
        return self.coordsList
        
    def asQgsPointList(self):
        """ Return the list of points.
        
        :rtype: list
        :returns: The points list.
        """
        return self.pointsList
        
    def asQgsFeature(self):
        """ Return the line feature. """
        lineFeature = qgis.core.QgsFeature()
        lineFeature.setGeometry(qgis.core.QgsGeometry.fromPolylineXY(self.pointsList))
        return lineFeature
        
    def asShapes(self):
        """ Return the line shape. """
        from .shapefileTools import Shapes
        lineFeature = self.asQgsFeature()
        lineShape = Shapes()
        lineShape.setShapeType(qgis.core.QgsWkbTypes.LineString)
        lineShape.setCoordRefSystem(self.coordRefSystem)
        lineShape.addFeature(lineFeature)
        return lineShape

class polygon(object):
    """ A class to represent polygon objects in the mesh. """
    def __init__(self, coordsList, coordRefSystem):
        """ Class initialisation method.
               
        :arg list coordsList: The list of coordinates which make up the points in the polygon.
        :arg coordRefSystem: The coordinate reference system to use. This can be a qgis.core.QgsCoordinateReferenceSystem object, or a string representing the CRS' name.
        """
        if type(coordRefSystem) != str and \
           type(coordRefSystem) != qgis.core.QgsCoordinateReferenceSystem:
            raise BadArguments('Error: Input agrument coordRefSystem is of'+\
                               ' incorrect type, expected a string or'+\
                               ' QgsCoordinateReferenceSystem.\n')
        elif type(coordRefSystem) == str:
            self.coordRefSystem_string = coordRefSystem
            self.coordRefSystem = qgis.core.QgsCoordinateReferenceSystem()
            self.coordRefSystem.createFromString(coordRefSystem)
            if not self.coordRefSystem.isValid():
                msg = 'Error: could not initialise coordinate'+\
                      ' reference system '+self.coordRefSystem_string+' .\n'
                raise BadArguments(msg)
        else:
            self.coordRefSystem = coordRefSystem
            self.coordRefSystem_string = coordRefSystem.authid()
        #TODO: give a helpful message for empty coordsList
        self.coordsList = coordsList
        # Construct a QgsPoint-list from the coordinate-pairs. While the
        # same information is present in self.coordsList, we also
        # generate the pointsList, as it naturally goes together
        # with the coordinate reference system variable; especially
        # when changing the coordinate reference system. 
        self.pointsList = []
        for coordPair in self.coordsList:
            point = qgis.core.QgsPoint()
            point.setX(coordPair[0])
            point.setY(coordPair[1])
            self.pointsList.append(point)
            
    def changeCoordRefSystem(self, targetCRS_string):
        """ Change the coordinate reference system.
        
        :arg str targetCRS_string: The name of the coordinate reference system to use.
        """
        raise Exception('Method implementation is incomplete.\n')
        
    def asCoordsList(self):
        """ Return the list of coordinates.
        
        :rtype: list
        :returns: The coordinates list.
        """
        return self.coordsList
        
    def asQgsPointList(self):
        """ Return the list of points.
        
        :rtype: list
        :returns: The points list.
        """
        return self.pointsList
        
    def asQgsPolygonFeature(self):
        """ Return the polygon feature. """
        polygonFeature = qgis.core.QgsFeature()
        polygonFeature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY([self.pointsList]))
        return polygonFeature
        
    def asPolygonShapes(self):
        """ Return the polygon shape. """
        from .shapefileTools import Shapes
        polygonFeature = self.asQgsPolygonFeature()
        polygonShape = Shapes()
        polygonShape.setShapeType(qgis.core.QgsWkbTypes.Polygon)
        polygonShape.setCoordRefSystemFromString(self.coordRefSystem_string)
        polygonShape.addFeature(polygonFeature)
        return polygonShape
        
    def writeShapefile(self, filename):
        """ Write the polygon shape file.
        
        :arg str filename: The desired name of the shape file.
        """
        self.asPolygonShapes().writeFile(filename)

class loxodromicLine(Line):
    """ Class of loxodromic lines: Lines at a constant bearing. """
    def __init__(self, startPoint, trueNorthBearing, numbPoints, loopArounds, endPoint, coordRefSystem_string):
        """ Initialise a new loxodromicLine object.
        
        :arg tuple startPoint: A pair (2-tuple) of floating point numbers, specifying
                    the coordinates of the start point on the
                    loxodrome. 
        :arg float trueNorthBearing: The angle, in degrees, between the
                    loxodrome and the great circle towards the North
                    pole. Positive in the clockwise direction.
        :arg int numbPoints: The number of points on loxodrome segment.
        :arg int loopArounds: The number the loxodrome 'goes around' the
                    globe.
        :arg tuple endPoint: A pair (2-tuple) of floating point numbers, specifying
                    the coordinates of the end point on the
                    loxodrome.
                    
        If the start and end points are specified, they are checked to
        ensure they actually lie on the loxodrome. However, the user
        is at liberty to specify only one coordinate of the end point,
        while assigning a numpy.nan to the other coordinate. This
        function will then calculate the missing cooridnate taking into
        account the loopAround argument. """
        
        #TODO: check point variable is of correct type.
        self.startPoint = startPoint
        startPointXi = self.startPoint[0]
        startPointEta = self.startPoint[1]
        self.trueNorthBearing = trueNorthBearing
        self.alpha = np.tan(np.radians(90. - trueNorthBearing))
        self.beta = startPointEta - self.alpha*startPointXi
        self.endPoint = endPoint
        endPointXi = self.endPoint[0]
        endPointEta = self.endPoint[1]
        #numbPoints must be greater than or equal to 2
        if numbPoints < 2:
            raise BadArguments('Cannot create a loxodromic line segment'+
                               ' with less than two points.')
        #Loop arounds must be a natural number.
        if type(loopArounds) != int and loopArounds < 0:
            raise BadArguments('The loop-around argument of the'+
                                  ' loxodrome must be a natural number')
        #Check that end point, if fully specified, actually lies on
        # the loxodrome.
        if not(np.isnan(endPointXi)) and\
           not(np.isnan(endPointEta)):
            if abs(endPointEta - \
                   self.alpha*(endPointXi + loopArounds*360.0) - \
                   self.beta) > 1e-12:
                raise BadArguments('End point appears not to lie on'+
                                   ' loxodrome.')
        #Calculate coordinates of end-point, if not fully specified.
        if np.isnan(endPointXi) and not(np.isnan(endPointEta)):
            endPointXi = ((endPointEta - self.beta)/self.alpha) -\
                          loopArounds*360.0
        if not(np.isnan(endPointXi)) and np.isnan(endPointEta):
            endPointEta = self.alpha*(endPointXi + loopArounds*360.0) +\
                          self.beta
        #Create points along equi-spaced intervals on loxodrome.
        coordsList = []
        delta_lon = ((endPointXi + loopArounds*360.0) - startPointXi)/\
                     (numbPoints - 1)
        #Delta longitude will be 0 if the loxodrome is drawn at a
        # true-North-bearing of 0,  so lets calculate Delta latitude
        # as well.
        delta_lat = (endPointEta - startPointEta)/(numbPoints - 1)
        if delta_lon != 0.0:
            for pointIndex in range(numbPoints):
                lon = startPointXi + pointIndex*delta_lon
                lat = self.alpha*lon + self.beta
                lon = np.sign(lon)*np.mod(abs(lon),360.0)
                #if Latitude is greater than 180 or smaller than -180,
                # the loxodrome has reached the pole and must stop.
                if lat >= 180.0:
                    coordsList.append([lon,180.0])
                if lat <=-180.0:
                    coordsList.append([lon,-180.0])
                else:
                    coordsList.append([lon,lat])
        else:
            for pointIndex in range(numbPoints):
                lon = startPointXi
                lat = startPointEta + pointIndex*delta_lat
                #if Latitude is greater than 180 or smaller than -180,
                # the loxodrome has reached the pole and must stop.
                if lat >= 180.0:
                    coordsList.append([lon,180.0])
                if lat <=-180.0:
                    coordsList.append([lon,-180.0])
                else:
                    coordsList.append([lon,lat])
        super(loxodromicLine, self).__init__(coordsList, coordRefSystem_string)

#class orthodromicLine(object): to be continued...

class blendedLoxodromes(Line):
    """ A class to represent blended loxodrome lines in the mesh. """
    def __init__(self, loxodrome1, loxodrome2, numbPoints):
        """ Class initialisation method. """
        #TODO: check both lines have same CRS
        coordRefSystem_string = loxodrome1.coordRefSystem_string
        loxodrome1CoordsList = loxodrome1.asCoordsList()
        loxodrome2CoordsList = loxodrome2.asCoordsList()
        blendingFactor = np.linspace(0.,1.,len(loxodrome1CoordsList))
        coordsList = []
        for pointCounter in range(numbPoints):
            xi = loxodrome1CoordsList[pointCounter][0]*(1.-blendingFactor[pointCounter]) + loxodrome2CoordsList[numbPoints-1-pointCounter][0]*blendingFactor[pointCounter]
            eta = loxodrome1CoordsList[pointCounter][1]*(1.-blendingFactor[pointCounter]) + loxodrome2CoordsList[numbPoints-1-pointCounter][1]*blendingFactor[pointCounter]
            coordsList.append([xi, eta])
        super(blendedLoxodromes, self).__init__(coordsList, coordRefSystem_string)

class Circle(Line):
    """ A class to represent circle objects in the mesh. """
    def __init__(self, centerPointXi, centerPointEta, radius, numbPoints, coordRefSystem_string):
        """ Class initialisation method. """
        #TODO: Describe the arguments.
        
        self.centerPointXi = centerPointXi
        self.centerPointEta = centerPointEta
        if radius < 0:
            raise badArguments('Circle radious must be positive.')
        else:
            self.radius = radius
        if numbPoints < 3:
            raise badArguments('Must allow at least three points on the circle periphery.')
        else:
            self.numbPoints = numbPoints
        # Generate the points around the circle. Starting from the East, going around in an anti-clockwise fashion.
        coordsList = []
        delta_alpha = 2.*np.pi/numbPoints
        for pointCounter in range(numbPoints-1):
            xi = self.centerPointXi + np.sin(pointCounter*delta_alpha)*self.radius
            eta = self.centerPointEta + np.cos(pointCounter*delta_alpha)*self.radius
            coordsList.append([xi, eta])
        coordsList.append(coordsList[0])
        super(Circle, self).__init__(coordsList, coordRefSystem_string)
        
    def asQgsPolygonFeature(self):
        """ Return the polygon feature. """
        polygonFeature = qgis.core.QgsFeature()
        polygonFeature.setGeometry(qgis.core.QgsGeometry.fromPolygonXY([self.pointsList]))
        return polygonFeature
        
    def asPolygonShapes(self):
        """ Return the polygon shape. """
        from .shapefileTools import Shapes
        polygonFeature = self.asQgsPolygonFeature()
        polygonShape = Shapes()
        polygonShape.setShapeType(qgis.core.QgsWkbTypes.Polygon)
        polygonShape.setCoordRefSystemFromString(self.coordRefSystem_string)
        polygonShape.addFeature(polygonFeature)
        return polygonShape

#class ellipse(object): to be continued...

class BadArguments(Exception):
    """ A class used to raise exceptions for bad arguments. """
    def __init__(self, message):
        self.message=message
    def __str__(self):
        return self.message
