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
import logging
import numpy
import datetime
from ..__init__ import __git_sha_key__
from ..__init__ import __version__
from GFD_basisChangeTools import GFD_basisChangeTools
LOG = logging.getLogger(__package__)

class Geometry(object):
    def __init__(self):
        """ Initialise a new Geometry object. """
        LOG.info("Initialising geometry...")
        self.pointIndex = 0
        self.lineIndex = 0
        self.lineLoopIndex = 0
        self.surfaceIndex = 0
        self.volumeIndex = 0
        self.points={}
        self.lines={}
        self.bsplines={}
        self.lineLoops={}
        self.surfaces={}
        self.volumes={}
        self.physicalLineIDs={}
        self.physicalSurfaceIDs={}
        self.isSpherical=False
        LOG.info("Done.")

    def addPoint(self, point):
        """ Add a single point to the geometry.
        
        :arg list point: A list of coordinates in 3D, e.g. [0,1,1]
        :rtype: int
        :returns: The index of the new point.
        """
        assert type(point) == list
        assert len(point) == 3
        self.points[self.pointIndex] = point
        self.pointIndex += 1
        return self.pointIndex-1

    def addPoints(self, points):
        """ Add multiple points to the geometry.
        Note that this simply loops over the points and calls addPoint on each of them.
        
        :arg list points: A list containing 3D coordinates specified as lists, e.g. [ [0,1,1], [2,1,1], [3,1,2] ]
        :rtype: None
        """
        for point in points:
            self.addPoint(point)
        return
        
    def removePoint(self, pointIndex):
        """ Remove a point from the points dictionary.
        
        :arg int pointIndex: The index of the point to remove.
        :rtype: None
        """
        assert pointIndex < (len(self.points))
        del self.points[pointIndex]

    def getPoint(self, pointIndex):
        """ Return a specific point.
        
        :arg int pointIndex: The index of the point in the points dictionary.
        :rtype: list
        :returns: The coordinates of the point in 3D space (in a list).
        """
        assert pointIndex < (len(self.points))
        return self.points[pointIndex]

    def addPointsFromDictionary(self, pointDictionary):
        for point in pointDictionary:
            self.addPoint(pointDictionary[point])

    def pointsCount(self):
        """ Count the number of points in the geometry.
        
        :rtype: int
        :returns: The number of points in the geometry.
        """
        return len(self.points)

    def addLineSegment(self, pointIndices):
        """ Add a new line segment.
        
        :arg list pointIndices: A list of indices corresponding to existing points in the geometry.
        :rtype: int
        :returns: The index of the new line segment.
        """
        assert type(pointIndices) == list
        for pointIndex in pointIndices:
            assert pointIndex < (len(self.points))
        self.lines[self.lineIndex] = pointIndices
        self.lineIndex += 1
        return self.lineIndex-1

    def removeLineSegment(self, lineIndex):
        """ Remove a specific line segment.
        
        :arg int lineIndex: The index of the line segment in the geometry.
        :rtype: None
        """
        assert lineIndex < (len(self.lines))
        del self.lines[lineIndex]

    def getLineSegment(self, lineIndex):
        """ Return an existing line segment in the geometry.
        
        :arg int lineIndex: The index of the line segment that is to be returned.
        :rtype: list
        :returns: The line segment, in the form of a list of points.
        """
        
        assert lineIndex < (len(self.lines))
        return self.lines[lineIndex]

    def addLineSegmentsFromDictionary(self, lineDictionary):
        """ Add one or more line segments to the geometry from a dictionary of lines.
        
        :arg dict lineDictionary: A dictionary containing one or more lines.
        :rtype: None
        :raises AssertionError: if the addition of a line segment fails.       
        """
        for line in lineDictionary:
            try:
                self.addLineSegment(lineDictionary[line])
            except AssertionError:
                raise AssertionError

    def lineSegmentsCount(self):
        """ Count the number of line segments in the geometry.
        
        :rtype: int
        :returns: The number of line segments in the geometry.
        """
        return len(self.lines)

    def addBspline(self, pointIndices):
        assert type(pointIndices) == list
        for pointIndex in pointIndices:
            assert pointIndex < (len(self.points))
        self.bsplines[self.lineIndex] = pointIndices
        self.lineIndex += 1
        return self.lineIndex-1

    def removeBspline(self, lineIndex):
        assert lineIndex < (len(self.bsplines))
        del self.bsplines[lineIndex]

    def getBspline(self, lineIndex):
        assert lineIndex < (len(self.bsplines))
        return self.bsplines[lineIndex]

    def bsplinesCount(self):
        return len(self.bsplines)

    def addBsplinesFromDictionary(self, bsplineDictionary):
        for bspline in bsplineDictionary:
            try:
                self.addBspline(bsplineDictionary[bspline])
            except AssertionError:
                raise AssertionError

    def addLineLoop(self, lineIndices):
        """ Add a line loop to the geometry.
        
        :arg list lineIndices: The indices of the lines making up the loop.
        :rtype: int
        :returns: The index of the new line loop.       
        """
        assert type(lineIndices) == list
        #TODO: check lines exist.
        self.lineLoops[self.lineIndex] = lineIndices
        self.lineIndex += 1
        return self.lineIndex-1

    def lineLoopsCount(self):
        """ Count the number of line loops in the geometry.
        
        :rtype: int
        :returns: The number of line loops in the geometry.
        """
        return len(self.lineLoops)

    def addLineLoopsFromDictionary(self, lineLoopsDictionary):
        """ Add one or more line loops to the geometry from a dictionary of line loops.
        
        :arg dict lineLoopsDictionary: A dictionary containing one or more line loops.
        :rtype: None
        """
        #TODO: check keys from lineLoops dictionary not already in
        # self.lineLoops.keys
        for lineLoop in list(lineLoopsDictionary.keys()):
            self.lineLoops[lineLoop] = lineLoopsDictionary[lineLoop]
            self.lineIndex += 1

    def addPlaneSurface(self, lineIndices):
        """ Add a plane surface to the geometry.
        
        :arg list lineIndices: The indices of the lines making up the plane surface. These should form a line loop.
        :rtype: int
        :returns: The index of the new plane surface.
        """
        assert type(lineIndices) == list
        #TODO: check if loops exist.
        self.surfaces[self.surfaceIndex] = lineIndices
        self.surfaceIndex += 1
        return self.surfaceIndex-1

    def surfacesCount(self):
        """ Count the number of surfaces in the geometry.
        
        :rtype: int
        :returns: The number of surfaces in the geometry.
        """
        return len(self.surfaces)

    def addSurfacesFromDictionary(self, surfaceDictionary):
        """ Add one or more plane surfaces to the geometry from a dictionary of surfaces.
        
        :arg dict surfaceDictionary: A dictionary containing one or more plane surfaces.
        :rtype: None
        """
        #TODO: check if keys in surfaceDictionary already in
        # self.surfaces
        for surface in surfaceDictionary:
            self.addPlaneSurface(surfaceDictionary[surface])

    def addPhysicalLineID(self, lineIndex, physicalID):
        """ Assign a physical ID to the line with index lineIndex.
        The specified line must already exist in the geometry,
        an exception is otherwise raised. 
        
        :arg int lineIndex: The index of the line to assign a physical ID to.
        :arg int physicalID: The physical ID to be assigned.
        :rtype: None
        """
        #TODO: check if line with lineIndex exists. Raise exception if it does not
        if self.physicalLineIDs.have_key(physicalID):
            self.physicalLineIDs[physicalID].append(lineIndex)
        else:
            self.physicalLineIDs[physicalID] = [lineIndex]

    def physicalLineIDsCount(self):
        """ Return the number of Physical line IDs.
        
        :rtype: int
        :returns: The number of Physical line IDs.
        """
        return len(self.physicalLineIDs)

    def addPhysicalLineIDsFromDictionary(self, physicalIDsDictionary):
        """Assign physical ID to lines in input dictionary. The
        dictionary-key is the physical-ID and the value is a list
        of the lines (IDs) to be assigned this physical ID. 
        
        :arg dict physicalIDsDictionary: The dictionary containing the lines to assign physical IDs to.
        :rtype: None
        """
        for physicalID in physicalIDsDictionary:
            lines = physicalIDsDictionary[physicalID]
            #TODO: check that the lines actually exist. Raise exception if it does not
            if physicalID in self.physicalLineIDs:
                self.physicalLineIDs[physicalID] += lines
            else:
                self.physicalLineIDs[physicalID] = lines

    def addPhysicalSurfaceID(self, surfaceIndex, physicalID):
        '''Assign a physical ID to the surface with lindex
        surfaceIndex. The specified surface must already exist
        in the geometry, an exception is otherwise raised. '''
        #TODO: check if surface with SurfaceIndex exists. Raise exception if it does not
        if self.physicalSurfaceIDs.have_key(physicalID):
            self.physicalSurfaceIDs[physicalID].append(surfaceIndex)
        else:
            self.physicalSurfaceIDs[physicalID] = [surfaceIndex]

    def physicalSurfaceIDsCount(self):
        """ Count the number of Physical surface IDs in the geometry.
        
        :rtype: int
        :returns: The number of Physical surface IDs in the geometry.
        """
        return len(self.physicalSurfaceIDs)

    def addPhysicalSurfaceIDsFromDictionary(self, physicalIDsDictionary):
        """ Assign physical ID to surfaces in input dictionary. The
        dictionary-key is the physical-ID and the value is a list
        of the surfaces (IDs) to be assigned this physical ID.
        
        :arg dict physicalIDsDictionary: The dictionary containing the surfaces to assign physical IDs to.
        :rtype: None
        """
        for physicalID in physicalIDsDictionary:
            surfaces = physicalIDsDictionary[physicalID]
            #TODO: check that the surfaces actually exist. Raise exception if it does not
            if physicalID in self.physicalSurfaceIDs:
                self.physicalSurfaceIDs[physicalID] += surfaces
            else:
                self.physicalSurfaceIDs[physicalID] = surfaces

    def addPolarSphere(self, radius=6.37101e+6):
        self.centerID  = self.addPoint([0,0,0])
        self.northPoleID = self.addPoint([0,0,radius])
        self.surfaceIndex += 1
        self.polarSphereID = self.surfaceIndex
        self.isSpherical=True
        return self.polarSphereID

    def writeGmshGeoFile(self, geoFilename):
        """ Write a Gmsh .geo file describing the geometry. 
        
        :arg str geoFileName: The desired name of the .geo file.
        """
        LOG.info('Writing file '+geoFilename)
        # Open file.
        gmshFile = open(geoFilename,'w')
        # Write origin and time-stamp
        gmshFile.write('// Gmsh .geo file produced by qmesh version '+__version__+ \
                       ' (git sha key: ' + str(__git_sha_key__)+').\n')
        gmshFile.write('// ' + datetime.datetime.now().strftime("Date,time file written (yyyy/mm/dd, hour:minute:second): %Y/%m/%d , %H:%M:%S") + '\n')
        # If on-the-sphere make the sphere
        if self.isSpherical:
            LOG.debug('Writing polar sphere definition to '+geoFilename)
            gmshFile.write('\n// Definition of sphere, for conversion of coords from polar-stereographic projection:\n')
            gmshFile.write('Point( '+str(self.centerID)+' ) = {' +
                                            str(self.points[self.centerID][0])+','+
                                            str(self.points[self.centerID][1])+','+
                                            str(self.points[self.centerID][2])+' };\n')
            gmshFile.write('Point( '+str(self.northPoleID)+' ) = {' +
                                            str(self.points[self.northPoleID][0])+','+
                                            str(self.points[self.northPoleID][1])+','+
                                            str(self.points[self.northPoleID][2])+' };\n')
            gmshFile.write('PolarSphere( '+str(self.polarSphereID)+' ) = {' +
                                            str(self.centerID)+','+
                                            str(self.northPoleID)+' };\n')
        # Write points. In case a polar-sphere has already been written, do not
        # write again the points used to define the sphere, and tranform the points
        # to the Arctic polar-stereographic projection.
        if self.pointsCount() > 0 and not(self.isSpherical):
            LOG.debug('Writing ' + str(self.pointsCount()) +
                    ' point definitions to ' + geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.pointsCount()) + ' points:\n')
            for index in list(self.points.keys()):
                gmshFile.write('Point( '+str(index)+' ) = {' +
                                            str(self.points[index][0])+','+
                                            str(self.points[index][1])+','+
                                            str(self.points[index][2])+' };\n')
        if self.pointsCount() > 3 and self.isSpherical:
            LOG.debug('Transforming point coordinates into Arctic Polar Stereographic CRS')
            LOG.debug('and writing ' + str(self.pointsCount()-2) +
                    ' point definitions to '+geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.pointsCount()-2) + ' points:\n')
            for index in list(self.points.keys()):
                if index == self.centerID or index == self.northPoleID:
                    pass
                else:
                    # Transform point coords from longitute-latitude-radius to polar-stereographic.
                    surface_radius = self.points[self.northPoleID][2]
                    lon_lat_rad_coords = [self.points[index][0], self.points[index][1], surface_radius]
                    polar_stereo_coords = GFD_basisChangeTools.lonlatradius_2_polarStereographic(
                        lon_lat_rad_coords, surface_radius)
                    # Write transformed coordinates to file.
                    gmshFile.write('Point( '+str(index)+' ) = {' +
                                            str(polar_stereo_coords[0])+','+
                                            str(polar_stereo_coords[1])+','+
                                            '0.0 };\n')
        # Write lines
        if self.lineSegmentsCount() > 0:
            LOG.debug('Writing ' + str(self.lineSegmentsCount()) +
                    ' line definitions to ' + geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.lineSegmentsCount()) +
                    ' line-segments:\n')
            for index in list(self.lines.keys()):
                pointsInLineString = str(self.lines[index])
                pointsInLineString = pointsInLineString.replace('[','{')
                pointsInLineString = pointsInLineString.replace(']','}')
                gmshFile.write('Line( '+str(index)+' ) = '+pointsInLineString+';\n')
        # Write bsplines
        if self.bsplinesCount() > 0:
            LOG.debug('Writing ' + str(self.bsplinesCount()) +
                    ' B-spline definitions to ' + geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.bsplinesCount()) +
                    ' B-splines:\n')
            for index in list(self.bsplines.keys()):
                pointsInBsplineString = str(self.bsplines[index])
                pointsInBsplineString = pointsInBsplineString.replace('[','{')
                pointsInBsplineString = pointsInBsplineString.replace(']','}')
                gmshFile.write('BSpline( '+str(index)+' ) = '+pointsInBsplineString+';\n')
        # Write line-loops
        if self.lineLoopsCount() > 0:
            LOG.debug('Writing ' + str(self.lineLoopsCount()) +
                    ' line-loop definitions to ' + geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.lineLoopsCount()) +
                    ' line-loops:\n')
            for index in list(self.lineLoops.keys()):
                entitiesInLoopString = str(self.lineLoops[index])
                entitiesInLoopString = entitiesInLoopString.replace('[','{')
                entitiesInLoopString = entitiesInLoopString.replace(']','}')
                gmshFile.write('Line Loop( '+str(index)+' ) = '+entitiesInLoopString+';\n')
        # Write surfaces
        if self.surfacesCount() > 0:
            LOG.debug('Writing ' + str(self.surfacesCount()) +
                    ' surface definitions to ' + geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.surfacesCount()) +
                    ' surfaces:\n')
            for index in list(self.surfaces.keys()):
                loopsString = str(self.surfaces[index])
                loopsString = loopsString.replace('[','{')
                loopsString = loopsString.replace(']','}')
                gmshFile.write('Plane Surface( '+str(index)+' ) = '+loopsString+';\n')
        # Write physical line IDs
        if self.physicalLineIDsCount() > 0:
            LOG.debug('Writing ' + str(self.physicalLineIDsCount()) +
                    ' physical line ID definitions to '+geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.physicalLineIDsCount()) +
                    ' physical line ID`s:\n')
            for index in list(self.physicalLineIDs.keys()):
                linesString = str(self.physicalLineIDs[index])
                linesString = linesString.replace('[','{')
                linesString = linesString.replace(']','}')
                gmshFile.write('Physical Line( '+str(index)+' ) = '+linesString+';\n')
        # Write Physical surface IDs
        if self.physicalSurfaceIDsCount() > 0:
            LOG.debug('Writing ' + str(self.physicalSurfaceIDsCount()) +
                    ' physical surface ID definitions to '+geoFilename)
            gmshFile.write('\n// Definitions of ' + str(self.physicalSurfaceIDsCount()) +
                    ' physical surface ID`s:\n')
            for index in list(self.physicalSurfaceIDs.keys()):
                surfacesString = str(self.physicalSurfaceIDs[index])
                surfacesString = surfacesString.replace('[','{')
                surfacesString = surfacesString.replace(']','}')
                gmshFile.write('Physical Surface( '+str(index)+' ) = '+surfacesString+';\n')
        # Done writing, close the file.
        gmshFile.close()
        LOG.info('Finished writing ' + geoFilename)

class Domain(object):
    """ Class encapsulating the Gmsh domain: geometry with
    mesh-size-metric fields. """
    def __init__(self):
        """ Initialise a new Domain object. """
        self.geometryLineShapes = None
        self.geometryPolygonShapes = None
        self.meshMetricField = None
        self.targetCoordRefSystem_string = None

    def setGeometry(self, geometryLineShapes, geometryPolygonShapes):
        """ Set the geometry of the domain. """
        self.geometryLineShapes = geometryLineShapes
        self.geometryPolygonShapes = geometryPolygonShapes
        #TODO: check both geometryLineShapes, geometryPolygonShapes have same isGlobal and southPoleCoordinates
        self.isGlobal = geometryLineShapes.isGlobal
        self.southPoleCoordinates = geometryLineShapes.southPoleCoordinates

    def setMeshMetricField(self, meshMetricRaster):
        ''' Set the mesh metric field. '''
        self.meshMetricField = meshMetricRaster

    def setTargetCoordRefSystem(self, targetCoordRefSystem_string, fldFillValue=None):
        ''' Set the target coordinate reference system (CRS).
        
        :arg str targetCoordRefSystem_string: The target CRS, specified as a string.
        :arg fldFillValue: The fill value for re-gridding purposes.
        :raises Exception: if the CRS argument is not a string.        
        '''
        # Check input argument is of type string
        if type(targetCoordRefSystem_string) != str:
            raise Exception('targetCoordRefSystem_string argument must be of type string')
        self.targetCoordRefSystem_string = targetCoordRefSystem_string
        self.fldFillValue = fldFillValue

    def getGeometry(self):
        """ Return the geometry object associated with the domain.
        
        :rtype: qmesh.Geometry
        :returns: The geometry object associated with the domain.        
        """
        return self.geometry

    def getMeshMetricField(self):
        """ Return the mesh metric field associated with the domain.
        
        :returns: The mesh metric field associated with the domain.        
        """
        return self.meshMetricField

    def writeGmshFiles(self, gmshGeoFileName, gmshFldFileName=None):
        """ Output the domain as a Gmsh file(s). """
        gmshGeometry = self.shp2geo()
        gmshGeometry.writeGmshGeoFile(gmshGeoFileName)
        if gmshFldFileName!=None:
            # Convert raster to target CRS and re-grid, if necessary.
            rasterCRS_string = self.meshMetricField.getCoordRefSystem().authid()
            if self.targetCoordRefSystem_string == 'PCC' and \
               rasterCRS_string == 'EPSG:4326':
                pass
            elif self.targetCoordRefSystem_string == rasterCRS_string:
                pass
            else:
                 self.meshMetricField.changeCoordRefSystem(self.targetCoordRefSystem_string)
                 xiMin, xiMax, etaMin, etaMax = self.meshMetricField.getRasterBounds()
                 numb_xiPoints, numb_etaPoints, delta_xi, delta_eta = self.meshMetricField.getRasterResolution()
                 if type(self.fldFillValue) == type(None):
                     raise Exception('The mesh-size raster must be re-gridded, no fill-value has been specified.')
                 self.meshMetricField.reGrid(xiMin, xiMax, etaMin, etaMax, \
                                             numb_xiPoints, numb_etaPoints, \
                                             fillValue=self.fldFillValue)
            # Write raster into appropriate file.
            self.meshMetricField.writefld(gmshFldFileName, self.targetCoordRefSystem_string)
            # Append field-defining statements to gmsh-geo-file.
            gmshGeoFile = open(gmshGeoFileName,'a')
            gmshGeoFile.write('\n//Mesh-metric definition:\n')
            gmshGeoFile.write(
               'Mesh.CharacteristicLengthExtendFromBoundary = 0;\n')
            gmshGeoFile.write('Field[0] = Structured;\n')
            gmshGeoFile.write('Field[0].FileName = "'+\
                          gmshFldFileName+'";\n')
            gmshGeoFile.write('Field[0].TextFormat = 1;\n')
            if self.targetCoordRefSystem_string == 'PCC':
                gmshGeoFile.write('Field[1] = LonLat;\n')
                gmshGeoFile.write('Field[1].IField = 0;\n')
            else:
                gmshGeoFile.write('Field[1] = MathEval;\n')
                gmshGeoFile.write('Field[1].F = "F0";\n')
            gmshGeoFile.write('Background Field = 1;\n')
            gmshGeoFile.close()

    def shp2geo(self, surfaceRadius=6.37101e+6):
        '''Convert Shapes-objects to gmsh-geometry objects.
        :arg float surfaceRadius: The radius of the sphere approximating Earth's surface
        :returns: Gmsh geometry object'''
        #Ensure features stored into line-shapes describe single-part lines.
        # Also, extract the 'PhysID' attribute, if it exists and append it
        # to appropriate list.
        lines = []
        linePhysIDList = []
        for feature in self.geometryLineShapes.featureList:
            featureGeometry = feature.geometry()
            shapeType = featureGeometry.wkbType()
            if not (shapeType == qgis.core.QgsWkbTypes.LineString or \
                    shapeType == qgis.core.QgsWkbTypes.LineString25D):
                msg = 'Method shp2geo can only be invoked on single-part line objects.\n'
                raise Exception(msg)
            lines.append(featureGeometry.asPolyline())
            featurePhysID = feature.attribute('PhysID')
            linePhysIDList.append(featurePhysID)
        #Ensure features stored into polygon-shapes describes a polygon.
        # Extract the 'PhysID' attribute, if it exists and append it
        # to appropriate list.
        polygons = []
        polygonPhysIDList = []
        for feature in self.geometryPolygonShapes.featureList:
            featureGeometry = feature.geometry()
            shapeType = featureGeometry.wkbType()
            if not (shapeType == qgis.core.QgsWkbTypes.Polygon or \
                    shapeType == qgis.core.QgsWkbTypes.Polygon25D):
                msg = 'Method shp2geo can only be invoked on single-part polygon objects.\n'
                raise Exception(msg)
            polygons.append(featureGeometry.asPolygon())
            featurePhysID = feature.attribute('PhysID')
            polygonPhysIDList.append(featurePhysID)
        # Transform geometries to desired coordinate reference system.
        # If target CRS is not PCC, we can use EPSG codes to transform to desired CRS
        if self.targetCoordRefSystem_string != 'PCC':
            self.geometryLineShapes.changeCoordRefSystem(self.targetCoordRefSystem_string)
            self.geometryPolygonShapes.changeCoordRefSystem(self.targetCoordRefSystem_string)
        # If target CRS is PCC, we must transform to EPSG:4326, in preparation
        # to the polar stereographic re-projection later when writing the mesh
        # generator geometry file (writeGmshGeoFile)
        elif self.targetCoordRefSystem_string == 'PCC':
            self.geometryLineShapes.changeCoordRefSystem('EPSG:4326')
            self.geometryPolygonShapes.changeCoordRefSystem('EPSG:4326')
        #Construct the following dictionaries (maps):
        # point dictionary: {pointID : point coordinates list (triplet; x,y,z)} (holds information for gmsh definition of points)
        #                   The points are extracted from the line definitions.
        # inverse point dictionary: {point coordinates string: pointID}
        # line ID-point ID dictionary: {lineID: list of pointIDs constituting the line} (holds information for gmsh definition of lines, bsplines etc.)
        # ring ID-point ID dictionary: {ringID: list of pointIDs constituting the ring}
        # point ID-ring IDs dictionary: {pointID: list of ringIDs where the key participates}
        # ring ID-line ID dictionary: {ringID: list of lineIDs constituting the ring} (holds information for gmsh definition of line-loops)
        # surface-point ID dictionary: {surfaceID: list of pointIDs constituting the surface}
        # surface-ring ID dictionary: {surfaceID: list of ringIDs constituting the surface, first ring is the outer boundary of the polygon}
        LOG.debug('Constructing dictionaries...')
        LOG.debug('        Constructing point dictionaries...')
        #Construct the point dictionary.
        # First construct a temporary point-dictionary to remove duplicate points
        temp_pointDictionary = {}
        pointsInLines = 0
        for line in lines:
            for point in line:
                pointsInLines += 1
                temp_pointDictionary[str([point.x(),point.y(),0.0])] = point
        LOG.debug('        Found '+str(pointsInLines)+' points - extracted from lines.')
        LOG.debug('        Found '+str(len(temp_pointDictionary))+' points after removing duplicates.')
        #Construct the point-dictionary now duplicates are gone.
        pointIndex = 0
        pointDictionary = {}
        for point in list(temp_pointDictionary.values()):
            pointDictionary[pointIndex] = [point.x(),point.y(), 0.0]
            pointIndex +=1
        #Construct the inverse-point-dictionary
        pointIndex = 0
        invPointDictionary = {}
        for point in list(pointDictionary.values()):
            invPointDictionary[str([point[0],point[1],0.0])] = pointIndex
            pointIndex +=1
        #Construct the line ID-point ID dictionary.
        LOG.debug('        Constructing lineID-pointID dictionary...')
        lineIndex = 0
        lineIDpointID_dictionary = {}
        for line in lines:
            pointIndices = []
            for point in line:
                pointIndex = invPointDictionary[str([point[0],point[1],0.0])]
                pointIndices.append(pointIndex)
            lineIDpointID_dictionary[lineIndex] = pointIndices
            lineIndex += 1
        #Construct the ring ID-point ID dictionary and the ring-ringID dictionary.
        LOG.debug('        Constructing ringID-pointID dictionary...')
        ringIndex = 0
        ringIDPointID_dictionary = {}
        ringID_dictionary = {}
        ringsInPolygonFile = 0
        for polygon in polygons:
            for ring in polygon:
                pointIndices = []
                for point in ring:
                    #Consistency check: see if point from polygon is indeed in the
                    # points maps (i.e. it exists in one of the lines). If not,
                    # The lines are not consistent with the polygons as they do not
                    # contain the same points.
                    #if str([point[0],point[1],0.0]) not in invPointDictionary.keys():
                    #    msg = 'Polygons and lines appear to be inconsistent, they do not'+\
                    #          ' contain the same points. Found point with coordinates '+\
                    #          str([point[0],point[1],0.0]) + ' in a polygon, but the same'+\
                    #          ' point does not exist in the given lines.\n'
                    #    raise Exception(msg)
                    pointIndex = invPointDictionary[str([point[0],point[1],0.0])]
                    pointIndices.append(pointIndex)
                ringID_dictionary[ringIndex] = ring
                ringIDPointID_dictionary[ringIndex] = pointIndices
                ringIndex += 1
                ringsInPolygonFile += 1
        LOG.debug('        Found '+str(len(polygons))+' polygons.')
        LOG.debug('        Found '+str(ringsInPolygonFile)+' ring(s) in polygons.')
        #Construct the pointID-ringIDs dictionary
        LOG.debug('        Constructing pointID-ringID dictionary...')
        pointIDRingID_dictionary = {}
        for ringIndex in list(ringIDPointID_dictionary.keys()):
            pointsOnRing = ringIDPointID_dictionary[ringIndex]
            for pointIndex in pointsOnRing:
                if pointIndex in list(pointIDRingID_dictionary.keys()):
                    pointIDRingID_dictionary[pointIndex].append(ringIndex)
                else:
                    pointIDRingID_dictionary[pointIndex] = [ringIndex]
        #Contruct ringID-lineID dictionary: Go through the lines -just
        # the first & last points- and find out which ring they participate
        # in. This assumes line and ring dictionaries contain the
        # same points!
        LOG.debug('        Constructing ringID-lineID dictionary...')

        ringIDLineID_dictionary = {}
        for ringKey in list(ringIDPointID_dictionary.keys()):
            lineList = []
            for lineKey in lineIDpointID_dictionary:
                if lineIDpointID_dictionary[lineKey][0] in ringIDPointID_dictionary[ringKey] and \
                   lineIDpointID_dictionary[lineKey][-1] in ringIDPointID_dictionary[ringKey]:
                    lineList.append(lineKey)
                ringIDLineID_dictionary[ringKey] = lineList
        #Make the ring-dictionary descibe the ring in a contiguous fashion.
        LOG.debug('        Making ringID-lineID dictionary describe rings in a contiguous fashion...')
        for ringKey in ringIDLineID_dictionary:
            lineList = ringIDLineID_dictionary[ringKey]
            #Pick first line in ring
            contiguousLinesList = [lineList.pop(0)]
            while len(lineList) > 0:
                #Get last line in contiguous line list
                currentLine = abs(contiguousLinesList[-1])
                if numpy.sign(contiguousLinesList[-1]) >= 0:
                    currentLineDirection = 1
                else:
                    currentLineDirection = -1
                #Get last point in the line.
                if currentLineDirection == 1:
                    currentPoint = lineIDpointID_dictionary[currentLine][-1]
                else:
                    currentPoint = lineIDpointID_dictionary[currentLine][0]
                #Locate the next line that contains that point,
                # and figure out if it is the first or the
                # last point in that line
                safetyCounter = len(lineList)
                for nextLine in lineList:
                    if currentPoint in lineIDpointID_dictionary[nextLine]:
                        if currentPoint == lineIDpointID_dictionary[nextLine][0] and currentLineDirection == 1:
                            contiguousLinesList.append(nextLine)
                        if currentPoint == lineIDpointID_dictionary[nextLine][0] and currentLineDirection == -1:
                            contiguousLinesList.append(nextLine)
                        if currentPoint == lineIDpointID_dictionary[nextLine][-1] and currentLineDirection == 1:
                            contiguousLinesList.append(-nextLine)
                        if currentPoint == lineIDpointID_dictionary[nextLine][-1] and currentLineDirection == -1:
                            contiguousLinesList.append(-nextLine)
                        lineList.remove(nextLine)
                        break
                    safetyCounter-=1
                    if safetyCounter == 0:
                        brokenPoint = pointDictionary[currentPoint]
                        msg = 'Error: Input geometry appears to contain'+\
                              ' open loops. Check geometry at point ('+\
                              str(brokenPoint[0])+', '+\
                              str(brokenPoint[1])+') (CRS '+\
                              self.targetCoordRefSystem_string+').'
                        raise BadGeometry(msg)
            ringIDLineID_dictionary[ringKey] = contiguousLinesList
        #Re-process dictionaries involving ring IDs: Gmsh wants loops (rings) to
        # be numbered with the same counter as other line types. However, so 
        # far we have numbered rings starting from 0. Create map oldIndex:newIndex
        # and use it to re-number the dictionaries.
        ringIDrenumbering = {}
        for ringID in list(ringIDLineID_dictionary.keys()):
            ringIDrenumbering[ringID] = ringID + len(lineIDpointID_dictionary)
        #Renumber entities in ring ID-point ID dictionary
        renumbered_ringIDPointID_dictionary = {}
        renumbered_ringID_dictionary = {}
        for ringID in list(ringIDPointID_dictionary.keys()):
            renumbered_ringIDPointID_dictionary[ringIDrenumbering[ringID]] = ringIDPointID_dictionary[ringID]
            renumbered_ringID_dictionary[ringIDrenumbering[ringID]] = ringID_dictionary[ringID]
        ringIDPointID_dictionary = renumbered_ringIDPointID_dictionary
        ringID_dictionary = renumbered_ringID_dictionary
        #Renumber entities in point ID-ring ID dictionary
        renumbered_pointIDRingID_dictionary = {}
        for pointID, ringIDs in zip(list(pointIDRingID_dictionary.keys()), list(pointIDRingID_dictionary.values())):
            renumbered_ringIDs=[]
            for old_ringID in ringIDs:
                renumbered_ringIDs.append(ringIDrenumbering[old_ringID])
            renumbered_pointIDRingID_dictionary[pointID] = renumbered_ringIDs
        pointIDRingID_dictionary = renumbered_pointIDRingID_dictionary
        #Renumber entities in ring ID-line ID dictionary
        renumbered_ringIDLineID_dictionary = {}
        for ringID in list(ringIDLineID_dictionary.keys()):
            renumbered_ringIDLineID_dictionary[ringIDrenumbering[ringID]] = ringIDLineID_dictionary[ringID]
        ringIDLineID_dictionary = renumbered_ringIDLineID_dictionary
        #Construct surfaceID-pointID dictionary.
        LOG.debug('        Constructing surfaceID-pointID dictionary...')
        #TBD?surfaceIDPointID_dictionary = {}
        surfaceIDRingID_dictionary = {}
        surfaceIndex = 0
        for polygon in polygons:
            #TBD?surfaceIDPointID_dictionary[surfaceIndex] = []
            surfaceIDRingID_dictionary[surfaceIndex] = []
            for ring in polygon:
                ringID = list(ringID_dictionary.keys())[list(ringID_dictionary.values()).index(ring)]
                surfaceIDRingID_dictionary[surfaceIndex].append(ringID)
                #TBD?pointsInRing = []
                #TBD?for point in ring:
                #TBD?    pointID = invPointDictionary[str([point[0],point[1],0.0])]
                #TBD?    pointsInRing.append(pointID)
                #TBD?surfaceIDPointID_dictionary[surfaceIndex].append(pointsInRing)
            surfaceIndex += 1
        #TBD?#Construct surfaceID-ringID dictionary.
        #TBD?if verbocityLevel > 2:
        #TBD?    sys.stdout.write('        Constructing surfaceID-ringID dictionary...\n')
        #TBD?    sys.stdout.flush()
        #TBD?surfaceIDRingID_dictionary = {}
        #TBD?for surfaceIndex in surfaceIDPointID_dictionary.keys():
        #TBD?    surfaceIDRingID_dictionary[surfaceIndex] = []
        #TBD?    for ring in surfaceIDPointID_dictionary[surfaceIndex]:
        #TBD?        ringID = ringID_dictionary.keys()[ringID_dictionary.values().index(ring)]
        #TBD?        surfaceIDRingID_dictionary[surfaceIndex].append(ringID)
        #Construct physical-line ID dictionary, if available
        LOG.debug('        Constructing Physical IDs dictionaries for lines...')
        physLineIDDictionary = {}
        if len(linePhysIDList) > 0:
            for line, lineIndex in zip(linePhysIDList, list(range(len(linePhysIDList)))):
                if str(line) != 'NULL' and str(line) != 'None':
                    if line in physLineIDDictionary:
                        physLineIDDictionary[line].append(lineIndex)
                    else:
                        physLineIDDictionary[line] = [lineIndex]
        #Construct physical-surface ID dictionary, if available
        LOG.debug('        Constructing Physical IDs dictionaries for surfaces...')
        physSurfaceIDDictionary = {}
        if len(polygonPhysIDList) > 0:
            for polygon, polygonIndex in zip(polygonPhysIDList, list(range(len(polygonPhysIDList)))):
                if polygon in physSurfaceIDDictionary:
                    physSurfaceIDDictionary[polygon].append(polygonIndex)
                else:
                    physSurfaceIDDictionary[polygon] = [polygonIndex]
        #Construct Geometry object from dictionaries
        LOG.debug('Constructing Geometry object...')
        gmshGeometry = Geometry()
        gmshGeometry.addPointsFromDictionary(pointDictionary)
        gmshGeometry.addBsplinesFromDictionary(lineIDpointID_dictionary)
        gmshGeometry.addLineLoopsFromDictionary(ringIDLineID_dictionary)
        gmshGeometry.addSurfacesFromDictionary(surfaceIDRingID_dictionary)
        gmshGeometry.addPhysicalLineIDsFromDictionary(physLineIDDictionary)
        gmshGeometry.addPhysicalSurfaceIDsFromDictionary(physSurfaceIDDictionary)
        #Add polar-sphere if output is in planet-centered cartesian system
        if self.targetCoordRefSystem_string == 'PCC':
            gmshGeometry.addPolarSphere(surfaceRadius)

        return gmshGeometry

    def gmsh(self, geoFilename=None, fldFilename=None, mshFilename=None, gmshAlgo='del2d', isMshFileBinary=False, verbosityLevel=None):
        '''Generate mesh using gmsh with specified options
        '''
        import subprocess
        from ..config import _subprocess_log_queue
        if geoFilename==None:
            geoFilename = 'domain.geo'
        if mshFilename==None:
            mshFilename = 'mesh.msh'
        # Write geo and fld files
        self.writeGmshFiles(geoFilename, fldFilename)
        # Compose gmsh terminal command
        gmshCommand = ['gmsh','-2','-algo',gmshAlgo]
        if verbosityLevel != None:
            try:
              if verbosityLevel not in [0,1,10,100]:
                msg = 'gmsh vebrosity level must be one of 0,1,10,100'
                raise Exception(msg)
            except:
                LOG.error(msg, exc_info=True)
                raise
            gmshCommand.append('-v')
            gmshCommand.append(str(verbosityLevel))
        if isMshFileBinary:
            gmshCommand.append('-bin')
        gmshCommand.extend([geoFilename,'-o',mshFilename])
        # Invoke gmsh command
        LOG.info('Meshing with gmsh')
        gmsh_hasErrors = False
        #Cannot use pipes to capture gmsh stdout, as it seems to buffer output
        # when a pipe is present. Using files instead 
        qmesh_pid = subprocess.os.getpid()
        gmsh_stdoutFileName = '/tmp/gmsh_stdout_'+str(qmesh_pid)
        gmsh_stdout = open(gmsh_stdoutFileName,'w')
        gmsh_stderrFileName = '/tmp/gmsh_stderr_'+str(qmesh_pid)
        gmsh_stderr = open(gmsh_stderrFileName,'w')
        gmsh_proc = subprocess.Popen(gmshCommand,
                                     stdout=gmsh_stdout,
                                     stderr=gmsh_stderr)
        gmsh_log_queue = _subprocess_log_queue(gmsh_proc, 'gmsh', gmsh_stdoutFileName, gmsh_stderrFileName)
        #Get a child logger (child of the qmesh logger) to output stdout and stderr from gmsh
        gmsh_logger = logging.getLogger('qmesh.mesh.gmsh')
        while True:
          if gmsh_log_queue.emptyQueue() and gmsh_proc.poll() is not None:
            break
          elif gmsh_log_queue.emptyQueue() and gmsh_proc.poll() is None:
            continue
          elif not gmsh_log_queue.emptyQueue():
            line = gmsh_log_queue.getLine()
            log_message = line.strip('\n').split()
            if log_message[0] == 'Info':
              gmsh_logger.info(' '.join(log_message[2:]))
            elif log_message[0] == 'Debug':
              gmsh_logger.debug(' '.join(log_message[2:]))
            elif log_message[0] == 'Warning':
              gmsh_logger.warning(' '.join(log_message[2:]))
            elif log_message[0] == 'Error':
              gmsh_logger.error(' '.join(log_message[2:]))
              gmsh_hasErrors = True
            else:
              gmsh_logger.info(' '.join(log_message))
        #Raise exception if gmsh failed
        try:
          if gmsh_proc.returncode != 0 or gmsh_hasErrors:
            msg = 'Encountered errors in producing mesh with gmsh command: '+' '.join(gmshCommand[0:])
            raise Exception(msg)
        except:
            LOG.error(msg, exc_info=True)
            raise

        gmsh_stdout.close()
        gmsh_stderr.close()
        # If the domain is global and the coordinate reference system
        # has been rotated, rotate the mesh in the opposite direction.
        if self.isGlobal and self.southPoleCoordinates != None:
            LOG.info('Rotating PCC coordinate Reference System.')
            rotateMsh(mshFilename, self.southPoleCoordinates, invertRotation=True)
        LOG.info('Mesh written to '+mshFilename)

class Mesh():
      def __init__(self):
          self.crs=None
          self.vertexCoordinates = []
          self.vertexLabels = []
          self.edgeLabels = []
          self.triangleEdges = []
          self.triangleLabels = []

      def readGmsh(self, gmshFileName, crs):
          '''Read a gmsh mesh, from a given file.
             targetCrs can be WTK, qgis.core.QgsCoordinateReferenceSystem()  or have the value 'PCC'
          '''
          #Set CRS to given value
          if crs=='PCC' or type(crs) == qgis.core.QgsCoordinateReferenceSystem:
              self.crs = crs
          else:
              coordinateReferenceSystem = qgis.core.QgsCoordinateReferenceSystem()
              coordinateReferenceSystem.createFromString(crs)
              self.crs = coordinateReferenceSystem
          #Read in mesh
          LOG.info("Reading mesh from "+gmshFileName)
          mshfile=open(gmshFileName, 'r')
          # Header section
          assert(mshfile.readline().strip()=="$MeshFormat")
          assert(mshfile.readline().strip()in["2 0 8", "2.1 0 8", "2.2 0 8"])
          assert(mshfile.readline().strip()=="$EndMeshFormat")
          #Read in vertices
          while mshfile.readline().strip() !="$Nodes":
              pass
          nodecount=int(mshfile.readline())
          if nodecount==0:
              raise Exception("ERROR: No nodes found in mesh.\n")
          for counter in range(nodecount): 
              # Node syntax
              line = mshfile.readline().split()
              # compare node id assigned by gmsh to consecutive node id (assumed by fluidity)
              if eval(line[0]) != counter+1:
                  raise Exception("ERROR: Nodes in gmsh .msh file must be numbered consecutively.\n")
              #Note, .msh files use 1-based numbering. 0-based numbering used here.
              self.vertexCoordinates.append(list(map(float,line[1:4])))
          assert(mshfile.readline().strip()=="$EndNodes")
          LOG.info('Read in '+str(nodecount)+' nodes.')
          #Read in edges and elements
          assert(mshfile.readline().strip()=="$Elements")
          elementcount=int(mshfile.readline())
          boundaryEdge_counter = 0
          quadrangle_counter = 0
          #TBD boundaryEdges=[]
          triangles=[]
          #TBD tets=[]
          #TBD quads=[]
          #TBD hexes=[]
          for counter in range(elementcount):
              element=mshfile.readline().split()
              if (element[1]=="1"):#Lines on boundaries of the domain
                  vertexIDs = [int(element[-2])-1 , int(element[-1])-1]
                  numberOfLabels = int(element[2])
                  elementLabels = list(map(int,element[3:3+numberOfLabels]))
                  boundaryEdge_counter += 1
                  #TBD boundaryEdges.append(map(int,element[-2:]+[element[3]]))
              elif (element[1]=="2"):
                  vertexIDs = [int(element[-3])-1 , int(element[-2])-1 , int(element[-1])-1]
                  numberOfLabels = int(element[2])
                  elementLabels = list(map(int,element[3:3+numberOfLabels]))
                  triangles.append(vertexIDs)
              elif(element[1]=="15"):
                  # Ignore point elements
                  pass
              else:
                  try:
                    raise Exception()
                  except Exception:
                    LOG.error("Unknown element type "+str(element[1])+
                              " in gmsh msh file "+gmshFileName, exc_info=True)
                    raise
          if boundaryEdge_counter > 0:
              LOG.info('Read in '+str(boundaryEdge_counter)+' 2-node edges on boundaries.')
          if len(triangles) > 0:
              LOG.info('Read in '+str(len(triangles))+' 3-node triangles.')
          if quadrangle_counter > 0:
              LOG.info('Read in '+str(len(quads))+' 4-node quadrangles')
          #TBU LOG.info('Read in '+str(len(tets))+' tets.')
          #TBU LOG.info('Read in '+str(len(hexes))+' hexes.')
          #The list "triangles" is a connectivity table, where for each element
          # the IDs of participating edges are listed. If we just draw all elemets
          # we will then draw each edge twice. Below we remove such duplications.
          # Create an edge:triangle dictionary first. Loop through the triangles
          # and accumulate edges in a separate list, where each edge is listed
          # only once. Note that below we assume each triangle edge is defined
          # from the two end-vertices, no "curved" triangles!
          LOG.info('Removing duplicate edge definitions...')
          tempTriangleEdges = []
          for triangle in triangles:
              tempTriangleEdges.append([triangle[0],triangle[1]])
              tempTriangleEdges.append([triangle[1],triangle[2]])
              tempTriangleEdges.append([triangle[2],triangle[0]])
          #Loop through edges list and sort each edge. Makes it easier to detect multiply defined edges
          #Operation list(set()) will remove duplicate entries, but cannot operate on list-of-lists, so
          #we must convert point ID lists into strings.
          triangleEdgeStrings = []
          for triangleEdge in tempTriangleEdges:
              triangleEdge.sort
              triangleEdgeStrings.append(str(triangleEdge))
          triangleEdgeStrings = list(set(triangleEdgeStrings))
          #triangleEdgeStrings is now a list specifying the edges of the mesh, in terms of point IDs.
          # But it is using string variables, so we now convert those into integers.
          for triangleEdgeString in triangleEdgeStrings:
              pointID_strings = triangleEdgeString[1:-1].split(',')
              triangleEdge = [int(pointID_strings[0]), int(pointID_strings[1]) ]
              self.triangleEdges.append(triangleEdge)
          LOG.info('Found '+str(len(triangleEdgeStrings))+' edges...')


      def writeGmsh(self, filename):
              raise Exception('NOT IMPLEMENTED')

      def writeShapefile(self, filename):
          from ..vector import shapefileTools
          if filename[-4:]!='.shp':
              filename += '.shp'
          LOG.info('Creating line-shapefile '+filename)
          #Create new shapefile object, loop trough triangle edges and add each
          # edge as a line.
          meshEdgeShapes = shapefileTools.Shapes()
          meshEdgeShapes.setShapeType(qgis.core.QgsWkbTypes.LineString)
          if type(self.crs) == str and self.crs != 'PCC':
              meshEdgeShapes.setCoordRefSystemFromString(self.crs)
          elif type(self.crs) == qgis.core.QgsCoordinateReferenceSystem:
              meshEdgeShapes.setCoordRefSystem(self.crs)
          elif self.crs == 'PCC':
              message = 'Cannot create shapefiles with "PCC" crs, '
              message += 'consider another projection for the mesh.'
              LOG.error(message)
              raise Exception(message)
          else:
              message = 'Mesh has incorrect Coordinate Reference System.'
              LOG.error(message)
              raise Exception(message)
          for edge in self.triangleEdges:
              pointsList = []
              for vertexCoordinates in [ self.vertexCoordinates[edge[0]],
                                         self.vertexCoordinates[edge[1]] ]:
                  point = qgis.core.QgsPoint()
                  point.setX(vertexCoordinates[0])
                  point.setY(vertexCoordinates[1])
                  pointsList.append(point)
              newFeature = qgis.core.QgsFeature()
              newFeature.setGeometry(qgis.core.QgsGeometry.fromPolyline(pointsList))
              meshEdgeShapes.addFeature(newFeature)
          #Write new features to shapefile.
          meshEdgeShapes.writeFile(filename)
          return

      def writeKml(self, filename):
          raise Exception('NOT IMPLEMENTED')

      def reProjectVertices(self, toCrs):
          '''Change the Coordinate Reference System used to define the vertex coordinates.
             targetCrs can be WTK, qgis.core.QgsCoordinateReferenceSystem()  or have the value 'PCC'
             It is highly recommended to avoid using this method and instead produce the mesh at the
             required CRS. However this method is provided for convenience.
          '''
          message = 'Avoid re-projecting mesh vertices. '
          message += 'Instead produce the mesh in the correct CRS/projection. '
          message += 'Re-projection can be used for visualisation purposes.'
          LOG.warning(message)
          #Set up the target Coordinate Reference System.
          if toCrs=='PCC' or type(toCrs) == qgis.core.QgsCoordinateReferenceSystem:
              targetCoordReferenceSystem = toCrs
          else:
              targetCoordReferenceSystem = qgis.core.QgsCoordinateReferenceSystem()
              targetCoordReferenceSystem.createFromString(toCrs)
          #Re-project vertex coordinates
          #If re-projecting onto PCC, then transform vertex coordinates from current projection
          # to authalic sphere, and then convert into PCC using GFD_basisChangeTools
          if targetCoordReferenceSystem == 'PCC' and type(self.crs) == qgis.core.QgsCoordinateReferenceSystem:
              lonLatCRS = qgis.core.QgsCoordinateReferenceSystem()
              lonLatCRS.createFromString('EPSG:4326')
              CoordTransformer = qgis.core.QgsCoordinateTransform(self.crs,
                                                                  lonLatCRS,
                                                                  qgis.core.QgsProject.instance())
              # Check if coordinate transformation is correctly initialised
              if not(CoordTransformer.isValid()):
                  message = 'Re-projection failed.'
                  message += ' Initialisation of coordinate'
                  message += ' transformation object failed.'
                  raise Exception(message)
              message = 'Converting mesh vertices from '+str(self.crs.description())
              message += ' to target reference coordinate system: PCC'
              LOG.info(message)
              for pointCounter in range(len(self.vertexCoordinates)):
                  sourceCRS_point = qgis.core.QgsPointXY(self.vertexCoordinates[pointCounter][0],
                                                       self.vertexCoordinates[pointCounter][1])
                  targetCRS_point = CoordTransformer.transform(sourceCRS_point)
                  pointCoordinates = GFD_basisChangeTools.lonlatradius_2_cartesian([targetCRS_point.x(), targetCRS_point.y(), 6371007.0])
                  self.vertexCoordinates[pointCounter] = pointCoordinates
                  #TBD self.vertexCoordinates[pointCounter] = [targetCRS_point.x(), targetCRS_point.y(), targetCRS_point.z()]
          #If re-projecting from PCC, then transform vertex coordinates from current projection
          # to authalic sphere using GFD_basisChangeTools, and then convert to target CRS
          elif type(targetCoordReferenceSystem) == qgis.core.QgsCoordinateReferenceSystem and self.crs=='PCC':
              lonLatCRS = qgis.core.QgsCoordinateReferenceSystem()
              lonLatCRS.createFromString('EPSG:4326')
              CoordTransformer = qgis.core.QgsCoordinateTransform(lonLatCRS,
                                                                  targetCoordReferenceSystem,
                                                                  qgis.core.QgsProject.instance())
              # Check if coordinate transformation is correctly initialised
              if not(CoordTransformer.isValid()):
                  message = 'Re-projection failed.'
                  message += ' Initialisation of coordinate'
                  message += ' transformation object failed.'
                  raise Exception(message)
              message = 'Converting mesh vertices from PCC'
              message += ' to target reference coordinate system:'
              message += str(targetCoordReferenceSystem.description())
              LOG.info(message)
              for pointCounter in range(len(self.vertexCoordinates)):
                  pointCoordinates = self.vertexCoordinates[pointCounter]
                  pointCoordinates = GFD_basisChangeTools.cartesian_2_lonlatradius(self.vertexCoordinates[pointCounter])
                  sourceCRS_point = qgis.core.QgsPointXY(pointCoordinates[0], pointCoordinates[1])
                  targetCRS_point = CoordTransformer.transform(sourceCRS_point)
                  self.vertexCoordinates[pointCounter] = [targetCRS_point.x(), targetCRS_point.y(), 0]
          elif type(targetCoordReferenceSystem) == qgis.core.QgsCoordinateReferenceSystem and \
               type(self.crs) == qgis.core.QgsCoordinateReferenceSystem:
              raise Exception('NOT IMPLEMENTED')
              if targetCoordReferenceSystem.toWkt != self.crs.toWkt:
                  LOG.info('Converting mesh vertices to target reference'+
                           ' coordinate system: '+str(targetCoordReferenceSystem.description()))
                  CoordTransformer = qgis.core.QgsCoordinateTransform(sourceCRS, 
                                                                      targetCoordReferenceSystem,
                                                                      qgis.core.QgsProject.instance())
                  # Check if coordinate transformation is correctly initialised
                  if not(CoordTransformer.isValid()):
                      raise Exception('Change of Coordinate Reference System'+\
                                  ' failed. Initialisation of coordinate'+\
                                  ' transformation object failed.\n')
              else:
                  LOG.info('Source and target CRS are the same. No projection performed')
          else:
              if targetCoordReferenceSystem==self.crs=='PCC':
                  LOG.info('Source and target CRS are the same("PCC"). No projection performed')
              else:
                  raise Exception('Could not perform projection: source and/or target coordinate reference systems might be badly defined.')
          #Re-assign correct CRS
          self.crs = targetCoordReferenceSystem


def rotateMsh(inputMshFileName, southPoleCoordinates, invertRotation=False, outputMshFileName=None):
    '''Function rotating the coordinate system a mesh is in. Only makes sence for meshes that are on a plante-centred-Cartesian system.'''
    import fileinput
    import os
    from ..vector import shapefileTools
    # Open output node file - we write it on the fly
    inputLineCounter = 0
    if outputMshFileName == None:
        output_file = inputMshFileName+'.tmp'
    else:
        output_file = outputMshFileName
    # write new node coordinates to new node file
    f_out = open(output_file,'w')
    foundEndNodes = False
    for inputLine in fileinput.input([inputMshFileName]):
        if inputLine.strip() == '$EndNodes':
            foundEndNodes = True
        if (inputLineCounter > 4) and not foundEndNodes:
            data = inputLine.split()
            x = float(data[1])
            y = float(data[2])
            z = float(data[3])
            rotated_x, rotated_y, rotated_z = shapefileTools.rotateCartesianBasis([x,y,z], southPoleCoordinates, invertRotation=invertRotation)
            f_out.write(data[0]+" ")
            f_out.write(str(rotated_x)+" ")
            f_out.write(str(rotated_y)+" ")
            f_out.write(str(rotated_z)+" \n")
        else:
            f_out.write(inputLine)
        inputLineCounter += 1
    f_out.close()
    #Tidy-up files
    if outputMshFileName == None:
        os.system('mv -f '+output_file+' '+inputMshFileName)
