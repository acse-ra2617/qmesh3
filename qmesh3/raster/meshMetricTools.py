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
import sys
import numpy as np

import logging
LOG = logging.getLogger(__package__)

class raster(object):
    """Class for storage and handling of raster data.

    The primary purpose of this class is to store the data
    from rater data-objects in a pythonic way: The raster
    point coordinates and data are stored as numpy arrays.
    Some other useful information about the raster is also
    stored: The bounding coordinates, the coordinate reference
    system, the number of points along each direction, and
    an indication if the raster is of global coverage or not. 

    explain raster grid is square in logical space, and structured.
    """
    def __init__(self):
        """ Initialise a new Raster object. Set all of the attributes to None. """
        self.xiMin = None
        self.xiMax = None
        self.etaMin = None
        self.etaMax = None
        self.numb_xiPoints = None
        self.numb_etaPoints = None
        self.coordRefSystem = None
        self.xiData = None
        self.etaData = None
        self.variableData = None
        self.isGlobal = False
        
    def setGridData(self, xiData, etaData, variableData):
        """ Set the local attributes with data for xi and eta, as well as variable data, and calculate the raster bounds.
        
        :arg xiData: The raster grid xi cooridinates.
        :arg etaData: The raster grid eta cooridinates.
        :arg variableData: The raster grid variable data.
        """
        self.xiData = xiData
        self.etaData = etaData
        self.variableData = variableData
        self.numb_xiPoints = xiData.shape[1]
        self.numb_etaPoints = etaData.shape[0]
        self.calculateRasterBounds()
        
    def getGridData(self):
        """ Return arrays containing variable and point coordinates
        of raster.

        This method returns the data of raster grid points as numpy arrays.
        Three arrays are returned, the first containing the first coordinate,
        the second containing the second coordinate and the third containing
        the data. All three arrays are numb_etaPoints x numb_xiPoints.

        :rtype: tuple of numpy arrays
        :returns: xiData: The raster grid xi cooridinates.
        :returns: etaData: The raster grid eta cooridinates.
        :returns: variableData: The raster grid variable data.
        """
        return self.xiData, self.etaData, self.variableData
        
    def getCoordinatesData(self):
        """ Return 1d-arrays containing point coordinates
        of raster.

        This method returns the coordinates of raster grid points as numpy
        one-dimensional arrays (vectors). Two arrays are returned, the first
        containing the first coordinate (length of numb_xiPoints) and the
        second containing the second coordinate (length of numb_etaPoints)

        :rtype: tuple of numpy arrays
        :returns: xiData: The raster grid xi cooridinates.
        :returns: etaData: The raster grid eta cooridinates.
        """
        return (self.xiData[0,:], self.etaData[:,0])
        
    def setCoordRefSystem(self, coordRefSystem):
        """ Set the coordinate reference system (CRS). This assumes that the CRS has
        not already been set.
        
        :arg qgis.core.QgsCoordinateReferenceSystem coordRefSystem: The coordinate reference system to use.
        :arg string coordRefSystem: The coordinate reference system to use.
        :rtype: None
        :raises: Exception if the coordinate reference system has already been set.
        """
        # Setting the CRS is not allowed if that property has
        # alredy been set.
        if self.coordRefSystem != None:

            msg = 'Setting the CRS is only allowed if that property'+\
                  ' has not been alredy set. Did you mean to use'+\
                  ' changeCoordRefSystem instead?\n'
            raise Exception(msg + str(coordRefSystem))
        # so user can supply a string, we can then create the crs
        if (isinstance(coordRefSystem, str)):
            input_string = coordRefSystem
            coordRefSystem = qgis.core.QgsCoordinateReferenceSystem()
            coordRefSystem.createFromString(input_string)


        # Check coordRefSystem is of type QgsCoordinateReferenceSystem
        # or None, then initialise variable.
        if coordRefSystem != None and  \
           type(coordRefSystem) != qgis.core.QgsCoordinateReferenceSystem:
            msg = 'Error: variable coordRefSystem is of incorrect type,'+\
                  str(type(coordRefSystem)+'\n')
            raise BadArguments(msg)
        self.coordRefSystem = coordRefSystem
        
    def changeCoordRefSystem(self, targetCRS_string):
        """ Change the coordinate reference system (CRS). Map the coordinates
        currently in the geometry onto the new CRS.
        
        :arg qgis.core.QgsCoordinateReferenceSystem targetCRS_string: The target coordinate reference system.
        :raises: Exception if the CRS has not already been initialised.
        :rtype: None
        """
        # Contruct Coordinate Reference Systems and transformer object.
        self.calculateRasterBounds()
        sourceCRS = self.getCoordRefSystem()
        targetCRS = qgis.core.QgsCoordinateReferenceSystem()
        if targetCRS_string != 'PCC':
            targetCRS.createFromString(targetCRS_string)
        else:
            targetCRS.createFromString('EPSG:4326')
        LOG.info('Converting raster to target reference'+
             ' coordinate system: '+str(targetCRS.description()))
        CoordTransformer = qgis.core.QgsCoordinateTransform(sourceCRS, targetCRS, qgis.core.QgsProject.instance())
        # Check if coordinate transformation is correctly initialised
        if not(CoordTransformer.isValid()):
            raise Exception('Change of Coordinate Reference System'+\
                        ' failed. Initialisation of coordinate'+\
                        ' transformation object failed.\n')
        # Loop through points, calculate their coordinates in target CRS
        for xiCounter in range(self.numb_xiPoints):
            for etaCounter in range(self.numb_etaPoints):
                xi = self.xiData[etaCounter, xiCounter]
                eta = self.etaData[etaCounter, xiCounter]
                sourceCRS_point = qgis.core.QgsPointXY(xi, eta)
                targetCRS_point = CoordTransformer.transform(sourceCRS_point)
                self.xiData[etaCounter, xiCounter] = targetCRS_point.x()
                self.etaData[etaCounter, xiCounter] = targetCRS_point.y()
        self.coordRefSystem = targetCRS
        # Update raster bounds
        self.calculateRasterBounds()
        return
        
    def getRasterPeripheryPoints(self):
        """ Build a list of the points in the periphery of the grid. Assumes raster grid is square in some CRS.
        
        :rtype: list
        :returns: The points in the periphery of the grid.
        """
        peripheryPoints=[]
        for xiCounter in range(self.numb_xiPoints):
            xi = self.xiData[0, xiCounter]
            eta = self.etaData[0, xiCounter]
            point = qgis.core.QgsPoint(xi, eta)
            peripheryPoints.append(point)
        for etaCounter in range(self.numb_etaPoints):
            xi = self.xiData[etaCounter, self.numb_xiPoints-1]
            eta = self.etaData[etaCounter, self.numb_xiPoints-1]
            point = qgis.core.QgsPoint(xi, eta)
            peripheryPoints.append(point)
        for xiCounter in reversed(list(range(self.numb_xiPoints))):
            xi = self.xiData[self.numb_etaPoints-1, xiCounter]
            eta = self.etaData[self.numb_etaPoints-1, xiCounter]
            point = qgis.core.QgsPoint(xi, eta)
            peripheryPoints.append(point)
        for etaCounter in reversed(list(range(self.numb_etaPoints))):
            xi = self.xiData[etaCounter, 0]
            eta = self.etaData[etaCounter, 0]
            point = qgis.core.QgsPoint(xi, eta)
            peripheryPoints.append(point)
        return peripheryPoints
        
    def calculateRasterBounds(self):
        """ Calculate the bounds of the raster object and store the results in xiMin, xiMax, etaMin and etaMax.
        
        :return: None
        """
        # Build a list of the points in the periphery of the grid
        peripheryPoints = self.getRasterPeripheryPoints()
        # Calculate raster bounds, go around the
        # list of transformed coordinates of raster periphery.
        point = peripheryPoints.pop(-1)
        xiMin = point.x()
        xiMax = point.x()
        etaMin = point.y()
        etaMax = point.y()
        for point in peripheryPoints:
            if point.x()>xiMax:
                xiMax=point.x()
            if point.x()<xiMin:
                xiMin=point.x()
            if point.y()>etaMax:
                etaMax=point.y()
            if point.y()<etaMin:
                etaMin=point.y()
        self.xiMin = xiMin
        self.xiMax = xiMax
        self.etaMin = etaMin
        self.etaMax = etaMax
        LOG.debug('Grid bounds: '+str(xiMin)+' to '+str(xiMax)+",")
        LOG.debug(str(etaMin)+' to '+str(etaMax))
        
    def getCoordRefSystem(self):
        """ Return the coordinate reference system. 
        
        :rtype: qgis.core.QgsCoordinateReferenceSystem
        :return: The coordinate reference system currently in use.
        """
        return self.coordRefSystem

    def getRasterBounds(self):
        """ Return the bounds of the raster object.
        
        :rtype: tuple
        :return: The bounds of the raster as a tuple: min(xi), max(xi), min(eta), and max(eta)
        """
        return self.xiMin, self.xiMax, self.etaMin, self.etaMax

    def getDataBounds(self):
        """ Return the bounds of the data values associated with the raster object.
        
        :rtype: tuple
        :return: The bounds of the data values as a 2-tuple: min(data values), max(data values)
        """
        if self.variableData == None:
            msg = 'ERROR: The bounds of the raster can only be calculated'+\
                  ' once the object has been populated with data. Object'+\
                  ' seems to be not initialised.\n'
            raise Exception(msg)
            return None, None
        else:
            dataMin = np.min(self.variableData)
            dataMax = np.max(self.variableData)
            return dataMin, dataMax

    def getVariableData(self):
        """ Return raster point-data as numpy array.

        This method returns the grid data as a numpy array, coordinates not 
        included.

        :rtype: A numpy array
        :returns: A numpy array containing the grid data of the raster. The array indices are arranged as [etaCoord, xiCoord].
        """
        return self.variableData

    def setVariableData(self, variableData):
        """Method assigning raster point-data from given numpy array.

        This method assigns the grid data from a given numpy array,
        coordinates not included.

        :arg numpy.array variableData: The data at the raster grid points
        :rtype: None
        :raises: Exception is raised if the shape of the input array does not match the shape of the raster.
        """
        if variableData.shape != (self.numb_etaPoints, self.numb_xiPoints):
            raise Exception('Error while assigning raster data: Given array'+
                            ' does not match raster grid.\n')
        else:
            self.variableData = variableData

    def getRasterResolution(self):
        ''' Calculate the resolution of the raster object.
        
        :rtype: tuple
        :return: A tuple containing: the number of xi points, the number of eta points, the xi resolution, and the eta resolution.
        '''
        # Calculate raster resolution
        delta_xi = float(abs(self.xiMax - self.xiMin)/(self.numb_xiPoints))
        delta_eta = float(abs(self.etaMax - self.etaMin)/(self.numb_etaPoints))
        return self.numb_xiPoints, self.numb_etaPoints, delta_xi, delta_eta

    def fromFile(self, inputFileName):
        """ Read in a raster from a file.
        
        :arg str inputFileName: The name of the raster file to read.
        """
        try:
            import PyQt5 as PyQt
        except ModuleNotFoundError:
            import PyQt4 as PyQt
        LOG.info('Reading raster file '+inputFileName)
        fileInfo = PyQt.QtCore.QFileInfo(inputFileName)
        if not fileInfo.exists():
            raise Exception('Error: File '+inputFileName+' does not exist.')
        baseName = fileInfo.baseName()
        rasterLayer = qgis.core.QgsRasterLayer(inputFileName, baseName)
        #self.rasterLayer = rasterLayer
        if not rasterLayer.isValid():
            raise Exception('Error: failed to load raster file '+inputFileName)
        # Extract coordinate reference system
        crs = rasterLayer.crs()
        self.setCoordRefSystem(crs)
        # Extract data and create numpy arrays
        layerRectangle = rasterLayer.extent()
        self.xiMin = layerRectangle.xMinimum()
        self.xiMax = layerRectangle.xMaximum()
        self.etaMin = layerRectangle.yMinimum()
        self.etaMax = layerRectangle.yMaximum()
        self.numb_xiPoints = rasterLayer.width()
        self.numb_etaPoints = rasterLayer.height()
        xiVector = np.linspace(self.xiMin, self.xiMax,
                               self.numb_xiPoints)
        etaVector = np.linspace(self.etaMin, self.etaMax,
                                self.numb_etaPoints)
        self.xiData, self.etaData = np.meshgrid(xiVector, etaVector)
        self.variableData = np.zeros(self.xiData.shape)
        for xiCounter in range(self.numb_xiPoints):
            for etaCounter in range(self.numb_etaPoints):
                xi = xiVector[xiCounter]
                eta = etaVector[etaCounter]
                gridPoint = qgis.core.QgsPointXY(xi,eta)
                sample = rasterLayer.dataProvider().identify(gridPoint, \
                       qgis.core.QgsRaster.IdentifyFormatValue)
                self.variableData[etaCounter, xiCounter] = \
                                   sample.results()[1]
        self.calculateRasterBounds()

    def fromNetCDF(self, inputFileName, crs_string=None, loadInVariableName=None):
        """ Read in a raster from a NetCDF file.
        
        :arg str inputFileName: The name of the NetCDF file to be read.
        :arg str crs_string: The name of the coordinate reference system. This is optional.
        :arg str loadInVariableName: The name of the variable in the input file to be read in. This is optional.
        """
        from scipy.io import netcdf
        LOG.info("Reading NetCDF file "+inputFileName)
        # Open input netCDF raster file.
        inputNetCDF = netcdf.NetCDFFile(inputFileName, 'r')
        # Read-in dimension names
        dimensionNames = list(inputNetCDF.dimensions.keys())
        # We expect exactly two dimensions. Raise exception otherwise
        if len(dimensionNames) != 2:
            msg='ERROR: File '+inputFileName+' does not have exactly two'+\
                ' dimensions.\n'
            raise Exception(msg)
        # Read-in variable names
        variableNames = list(inputNetCDF.variables.keys())
        # Get rid of the "dimension-variables" from the variableNames list
        for dimensionName in dimensionNames:
            variableNames.remove(dimensionName)
        # If user has specified a named variable to load in, check it actually
        # is in the list. Otherwise load the first varible found in file.
        if loadInVariableName!=None and loadInVariableName not in variableNames:
            msg='ERROR: Variable named '+loadInVariable+\
                ' not found in file '+inputFileName+'\n'
            raise Exception(msg)
        else:
            loadInVariableName=variableNames[0]
        # Load-in the variable
        variable = inputNetCDF.variables[loadInVariableName]
        self.variableData = variable[:]
        # Find out the order the dimensions are defined in
        variableShape = self.variableData.shape
        if variableShape[0] == inputNetCDF.variables[dimensionNames[0]].shape:
            xiID = 0
            etaID = 1
        else:
            xiID = 1
            etaID = 0
        # Load in the dimensions
        xiDimensionName = dimensionNames[xiID]
        xiVariable = inputNetCDF.variables[xiDimensionName]
        xiVector = xiVariable[:]
        etaDimensionName = dimensionNames[etaID]
        etaVariable = inputNetCDF.variables[etaDimensionName]
        etaVector = etaVariable[:]
        self.xiData, self.etaData = np.meshgrid(xiVector, etaVector)
        self.xiMin = min(xiVector)
        self.xiMax = max(xiVector)
        self.etaMin = min(etaVector)
        self.etaMax = max(etaVector)
        self.numb_xiPoints = len(xiVector)
        self.numb_etaPoints = len(etaVector)
        self.calculateRasterBounds()
        # Look for a crs variable. If it does not exist, set coordinate reference system
        # from input argument of this method. If the user did not specify a coordinate
        # reference system, then default to EPSG:4326
        if 'crs' in variableNames:
            crsVariable = inputNetCDF.variables['crs']
            crsWellKnownText = getattr(crsVariable, 'spatial_ref')
            crs = qgis.core.QgsCoordinateReferenceSystem()
            crs.createFromWkt(crsWellKnownText)
            # If user has specified a crs while using this method
            # and it does not match that found in the file, inform
            # them of the mis-match and default to that found in file
            if crs_string!=None:
                if crs.authid()!=crs_string:
                    msg = '    WARNING: While reading file '+inputFileName+\
                          ' the following Coordinate-reference-system'+\
                          ' specification was found: '+crsWellKnownText+\
                          ' However, it does not match the given'+\
                          ' specification ('+crs_string+'). Defaulting to'+\
                          ' specification found in-file.'
                    LOG.warning(msg)
        elif 'crs' not in variableNames and crs_string!=None:
            crs = qgis.core.QgsCoordinateReferenceSystem()
            crs.createFromString(crs_string)
        else:
            #TBDcrs = qgis.core.QgsCoordinateReferenceSystem()
            #TBDcrs.createFromString('EPSG:4326')
            msg = 'Cannot find a Coordinate Reference System associated'+\
                  ' with raster file '+inputFileName+' . Please specify a'+\
                  ' CRS explicitly.\n'
            raise Exception(msg)
        self.setCoordRefSystem(crs)
        inputNetCDF.close()

    def fromIGF(self, inputFileName, crs_string=None, NODATA_value=np.nan):
        """ Read a raster from an IGF file. 
        
        :arg str inputFileName: The name of the IGF file to be read.
        :arg str crs_string: The name of the coordinate reference system. This is optional.
        :arg NODATA_value: The value that should replace NaNs. This is optional.
        """
        # Open file and extract grid specification data as strings
        inputFile = open(inputFileName, 'r')
        gridSpecification = {'ncols':0, 'nrows':0,
                             'xllcorner':0, 'yllcorner':0,
                             'cellsize':0, 'NODATA_value':np.nan}
        for line in inputFile:
            if line.split()[0] in list(gridSpecification.keys()):
                gridSpecification[line.split()[0]]=line.split()[1]
        inputFile.close()
        # Convert grid specification data into appropriate types
        gridSpecification['ncols'] = int(gridSpecification['ncols'])
        if gridSpecification['ncols'] <= 0:
            raise Exception('Error: it appears file '+inputFileName+\
                            ' describes a grid with '+str(gridSpecification['ncols'])+\
                            ' columns. Cannot proceed with constructing raster.\n')
        gridSpecification['nrows'] = int(gridSpecification['nrows'])
        if gridSpecification['nrows'] <= 0:
            raise Exception('Error: it appears file '+inputFileName+\
                            ' describes a grid with '+str(gridSpecification['nrows'])+\
                            ' rows. Cannot proceed with constructing raster.\n')
        gridSpecification['xllcorner'] = float(gridSpecification['xllcorner'])
        gridSpecification['yllcorner'] = float(gridSpecification['yllcorner'])
        gridSpecification['cellsize'] = float(gridSpecification['cellsize'])
        if gridSpecification['cellsize'] <= 0:
            raise Exception('Error: it appears file '+inputFileName+\
                            ' describes a grid with '+str(gridSpecification['cellsize'])+\
                            ' sized-cells. Cannot proceed with constructing raster.\n')
        gridSpecification['NODATA_value'] = float(gridSpecification['NODATA_value'])
        # Read-in variable data
        variableData = np.zeros((gridSpecification['nrows'], gridSpecification['ncols']))
        variableData = np.loadtxt(inputFileName, skiprows=6)
        variableData = np.flipud(variableData)
        # Reconstruct the raster data
        xiVector = np.linspace(
                gridSpecification['xllcorner'],
                gridSpecification['cellsize']*gridSpecification['ncols']+gridSpecification['xllcorner'],
                gridSpecification['ncols'])
        self.xiMin = min(xiVector)
        self.xiMax = max(xiVector)
        self.numb_xiPoints = len(xiVector)
        etaVector = np.linspace(
                gridSpecification['yllcorner'],
                gridSpecification['cellsize']*gridSpecification['nrows']+gridSpecification['yllcorner'],
                gridSpecification['nrows'])
        self.etaMin = min(etaVector)
        self.etaMax = max(etaVector)
        self.numb_etaPoints = len(etaVector)
        self.xiData, self.etaData = np.meshgrid(xiVector, etaVector)
        # Assign values to nans
        nanRows, nanColumns = np.where(variableData == gridSpecification['NODATA_value'])
        for rindex, cindex in zip(nanRows, nanColumns):
            variableData[rindex, cindex] = NODATA_value
        self.variableData = variableData
        self.calculateRasterBounds()
        # Set the coordinate reference system of the raster.
        if crs_string != None:
            crs = qgis.core.QgsCoordinateReferenceSystem()
            crs.createFromString(crs_string)
            self.setCoordRefSystem(crs)
        else:
            msg = 'Cannot find a Coordinate Reference System associated'+\
                  ' with raster file '+inputFileName+' . Please specify a'+\
                  ' CRS explicitly.\n'
            raise Exception(msg)

    def writeNetCDF(self, outputFileName,
                    xiName='x', xiAttributes={'long_name':'x coordinate'}, xiType=np.float64,
                    etaName='y', etaAttributes={'long_name':'y coordinate'}, etaType=np.float64,
                    variableName='z', variableAttributes={'long_name':'z variable'}, variableType=np.float64,
                    globalAttributes={'title':'z(x,y)'}):
        """ Write the raster to a NetCDF file. """
        from scipy.io import netcdf
        LOG.info('Writing NetCDF file '+outputFileName)
        #Open output netCDF raster file.
        outputNetCDF = netcdf.NetCDFFile(outputFileName, 'w')
        #Create the dimensions.
        outputNetCDF.createDimension(xiName, self.numb_xiPoints)
        outputNetCDF.createDimension(etaName, self.numb_etaPoints)
        #Output 'dimension variables' and their attributes,
        xiData, etaData = self.getCoordinatesData()
        if xiType == np.float64:
            xiVariable = outputNetCDF.createVariable(xiName, 'd', (xiName,))
            xiVariable[:] = np.float64(xiData)
        elif xiType == np.float32:
            xiVariable = outputNetCDF.createVariable(xiName, 'f', (xiName,))
            xiVariable[:] = np.float32(xiData)
        else:
            raise Exception('Output of xi-variable as given data-type to NetCDF'+\
                            ' file is not supported. Please choose amongst'+\
                            ' numpy.float64 and numpy.float32 .\n')
        for attribute in xiAttributes:
            setattr(xiVariable, attribute, xiAttributes[attribute])
        if etaType == np.float64:
            etaVariable = outputNetCDF.createVariable(etaName, 'd', (etaName,))
            etaVariable[:] = np.float64(etaData)
        elif etaType == np.float32:
            etaVariable = outputNetCDF.createVariable(etaName, 'f', (etaName,))
            etaVariable[:] = np.float32(etaData)
        else:
            raise Exception('Output of eta-variable as given data-type to NetCDF'+\
                            ' file is not supported. Please choose amongst'+\
                            ' numpy.float64 and numpy.float32 .\n')
        for attribute in etaAttributes:
            setattr(etaVariable, attribute, etaAttributes[attribute])
        #Output the raster variable and its attributes.
        if variableType == np.float64:
            variable = outputNetCDF.createVariable(variableName, 'd', (etaName, xiName))
            variable[:] = np.float64(self.variableData)
        elif variableType == np.float32:
            variable = outputNetCDF.createVariable(variableName, 'f', (etaName, xiName))
            variable[:] = np.float32(self.variableData)
        else:
            raise Exception('Output of variable as given data-type to NetCDF'+\
                            ' file is not supported. Please choose amongst'+\
                            ' numpy.float64 and numpy.float32 .\n')
        for attribute in variableAttributes:
            setattr(variable, attribute, variableAttributes[attribute])
        #Finally output the CRS variable and its attributes.
        if self.getCoordRefSystem() != None : 
            crsVariable = outputNetCDF.createVariable('crs', 'c', ())
            setattr(crsVariable, 'spatial_ref', str(self.getCoordRefSystem().toWkt()))
        #Set the global attributes
        for attribute in globalAttributes:
            setattr(outputNetCDF, attribute, globalAttributes[attribute])
        outputNetCDF.flush()
        outputNetCDF.close()

    def writeFluidityLonLatBathy(self, outputFileName,
                                 outputGridLon=None, outputGridLat=None):
       """ Write raster in a NetCDF file, formatted as Fluidity bathymetry.

       A Fluidity bathymetry NetCDF file is a structured grid in lon-lat
       (geographic coordinates). Furthermore the output NetCDF file is
       formated in a way that it can be read-in by Fluidity. The coordinate
       reference system of the raster is here checked, and if the grid is not
       in geographic coordinates it is reprojected, in which case the user
       must also specify the output grid coordinates arguments (see below).

       :arg str outputNetCDF_filename: A string storing the filename (with extension) of the output file.
       """
       xiName = 'lon'
       xiAttributes = {'long_name': 'longitude',
                       'units': 'degrees_east',
                       'actual_range': (self.xiMin, self.xiMax)
                      }
       etaName = 'lat'
       etaAttributes = {'long_name': 'latitude',
                        'units': 'degrees_north',
                        'actual_range': (self.etaMin, self.etaMax)
                       }
       variableName = 'z'
       variableAttributes = {
                        'long_name': 'z',
                        'actual_range': (np.nanmin(self.variableData), np.nanmax(self.variableData)),
                        '_FillValue': 'NaNf'
                       } 
       globalAttributes = {
                        'Conventions': 'COARDS/CF-1.0'
                       }
       self.writeNetCDF(outputFileName,
                    xiName, xiAttributes, xiType=np.float64,
                    etaName=etaName, etaAttributes=etaAttributes, etaType=np.float64,
                    variableName=variableName, variableAttributes=variableAttributes, variableType=np.float32,
                    globalAttributes=globalAttributes)

    def writeInfo(self):
        """ Print out information about the CRS, and the xi and eta points.
        
        :return: None
        """
        crs = self.getCoordRefSystem()
        if crs.description()!='' :
            LOG.debug('Found '+crs.description()+' ('+str(crs.authid())+') coordinate reference system.')
        elif crs.toWkt()!='' :
            LOG.debug('Found coordinate reference system with WKT: '+crs.toWkt())
        else:
            LOG.debug('Found no coordinate reference system.')
        if crs.description()!='':
            LOG.debug('Coordinate reference system has WKT: '+crs.toWkt())
        LOG.debug('Grid extents:')
        LOG.debug('xi: '+str(self.numb_xiPoints)+' points across '+str(self.xiMin)+' to '+str(self.xiMax))
        LOG.debug('eta: '+str(self.numb_etaPoints)+' points across '+str(self.etaMin)+' to '+str(self.etaMax))
        
    def writepng(self, outputFileName):
        """ Write the raster file to a Portable Network Graphics (PNG) image file.
        
        :arg str outputFileName: The desired filename of the PNG file.
        :return: None
        """
        import matplotlib.pyplot as plt
        plot = plt.pcolormesh(self.xiData, self.etaData, self.variableData)
        plt.colorbar()
        plt.title(self.coordRefSystem.authid())
        plt.savefig(outputFileName, format='png')
        plt.close()
        
    def writefld(self, outputFileName, targetCRS_string):
        """ Function for writing a gmsh field file, from the given
        raster-layer.

        :arg str outputFileName: The desired name of the gmsh file to be written.
        :arg str targetCRS_string: The name of the target coordinate reference system.
        :return: None
        """
        #TODO: Check grid is aligned with grid-coordinate-directions.
        #TODO: Check check target CRS is consistent with self.CRS
        # Open output file and write header
        LOG.info('Writing file '+outputFileName)
        numb_xiPoints, numb_etaPoints, delta_xi , delta_eta = self.getRasterResolution()
        outputFile=open(outputFileName,'w')
        if targetCRS_string == 'PCC':
            outputFile.write(str(np.radians(self.xiData[0,0]))+" "+\
                         str(np.radians(self.etaData[0,0]))+" 0\n")
            outputFile.write(
              str(np.radians(delta_xi))+" "+ str(np.radians(delta_eta))+" 1\n")
        else:
            outputFile.write(str(self.xiData[0,0])+" "+str(self.etaData[0,0])+" 0\n")
            outputFile.write(str(delta_xi)+" "+ str(delta_eta)+" 1\n")
        outputFile.write(str(numb_xiPoints)+" "+str(numb_etaPoints)+" 1\n")
        
        # Write data to output file, then close.
        for xiCounter in range(numb_xiPoints):
            for etaCounter in range(numb_etaPoints):
                # output 4dp floats with nans replaced with 0.000 to be compatible with GMSH
                outputFile.write(("%.4f\n" % self.variableData[etaCounter, xiCounter]).replace("nan", "0.000"))
        outputFile.close()
        LOG.info('Done writing file '+outputFileName)
        
    def getRasterPartitionPeripheryPoints(self, low_xiIndex, high_xiIndex, \
                                                low_etaIndex, high_etaIndex):
        """ Get the partition periphery points. """
        partitionPeripheryPoints=[]
        for partition_xiCounter in range(low_xiIndex, high_xiIndex+1):
            partitionPeripheryPoints.append( qgis.core.QgsPoint(
                             self.xiData[low_etaIndex, partition_xiCounter],
                             self.etaData[low_etaIndex, partition_xiCounter]))
        for partition_etaCounter in range(low_etaIndex, high_etaIndex+1):
            partitionPeripheryPoints.append( qgis.core.QgsPoint(
                             self.xiData[partition_etaCounter, high_xiIndex],
                             self.etaData[partition_etaCounter, high_xiIndex]))
        for partition_xiCounter in reversed(list(range(low_xiIndex, high_xiIndex+1))):
            partitionPeripheryPoints.append( qgis.core.QgsPoint(
                             self.xiData[high_etaIndex, partition_xiCounter],
                             self.etaData[high_etaIndex, partition_xiCounter]))
        for partition_etaCounter in reversed(list(range(low_etaIndex, high_etaIndex+1))):
            partitionPeripheryPoints.append( qgis.core.QgsPoint(
                             self.xiData[partition_etaCounter, low_xiIndex],
                             self.etaData[partition_etaCounter, low_xiIndex]))
        return partitionPeripheryPoints

    def reGrid(self, new_xiMin, new_xiMax, new_etaMin, new_etaMax, \
               new_numb_xiPoints, new_numb_etaPoints, \
               fillValue=0.0):
        """ Re-grid the mesh metric from the old grid to the new regular grid aligned with the CRS.
        
        :arg new_xiMin: The minimum xi coordinate.
        :arg new_xiMax: The maximum xi coordinate.
        :arg new_etaMin: The minimum eta coordinate.
        :arg new_etaMax: The maximum eta coordinate.
        :arg new_numb_xiPoints: The  number of xi points.
        :arg new_numb_etaPoints: The number of eta points.
        """
        from scipy.interpolate import griddata

        # Calculate point-coordinates of new grid
        new_xiList = np.linspace(new_xiMin, new_xiMax, new_numb_xiPoints)
        new_etaList = np.linspace(new_etaMin, new_etaMax, new_numb_etaPoints)
        newxiData, newetaData = np.meshgrid(new_xiList, new_etaList)

        oldM=self.xiData.shape[0]
        oldN=self.etaData.shape[1]
        newM=new_numb_xiPoints
        newN=new_numb_etaPoints
        xrange1=new_xiMin
        xrange2=new_xiMax
        yrange1=new_etaMin
        yrange2=new_etaMax
        msg='Regridding mesh metric from '+str(oldM)+'x'+str(oldN)+' grid to '+\
            'new grid aligned with the CRS, size '+str(newM)+'x'+str(newN)+', with'+\
            'spatial extent ['+str(xrange1)+','+str(xrange2)+']x['+str(yrange1)+','+str(yrange2)+']'
            
        LOG.info(msg)

        self.variableData = griddata(np.vstack((self.xiData.flatten(), self.etaData.flatten())).T,
            self.variableData.flatten(), 
            np.vstack((newxiData.flatten(), newetaData.flatten())).T,
            method = 'linear',
            fill_value = fillValue).reshape(newxiData.shape)
        self.xiData = newxiData
        self.etaData = newetaData

        self.numb_xiPoints = new_numb_xiPoints
        self.numb_etaPoints = new_numb_etaPoints
        self.calculateRasterBounds()

    def setGlobal(self, isGlobal=True):
        """ Set the boolean value of isGlobal.

        isGlobal helps to enable the correct differentiation of global
        bathymetries (+/- 180 longitude line).
        
        :arg bool isGlobal: The value that self.isGlobal should take.
        """
        self.isGlobal = isGlobal

    def deepcopy(self):
        """ Make a copy of the raster object.
        
        :rtype: raster
        :returns: A copy of the current raster object.
        """
        rasterCopy = raster()
        xiData, etaData, variableData = self.getGridData()
        rasterCopy.setGridData(xiData, etaData, variableData)
        rasterCopy.setCoordRefSystem(self.getCoordRefSystem())
        return rasterCopy

    def gradientMagnitude(self, targetCRS_string=None,
                          surfaceRadius=None):
        """ Calculate the raster gradient magnitude. 
        
        :arg str targetCRS_string: An optional argument to specify the coordinate reference
                                   system on which the magnitude is calculated.
        :arg float surfaceRadius: The radius of the surface.
        """
        # If the user specified a different CRS, get the data in that CRS
        inputRaster = self.deepcopy()
        xiData_degrees, etaData_degrees, variableData = self.getGridData()
        numb_xiPoints, numb_etaPoints, delta_xi, delta_eta = self.getRasterResolution()
        # When the output CRS is EPSG:4326 a 'strict' calculation of the gradient 
        # is done: The gradient is a vector quantity and as such it complies to
        # tensor definition. Otherwise the gradient is calculated from the partial
        # derivatives w.r.t. the coordinate directions alone.
        if targetCRS_string == None:
            targetCRS_string=self.coordRefSystem.authid()
        if targetCRS_string == 'EPSG:4326':
            # Make sure input CRS is lon-lat
            if inputRaster.getCoordRefSystem().authid() != 'EPSG:4326':
                inputRaster.changeCoordRefSystem('EPSG:4326')
            # Input coordinates might be in degrees, but here needed in radians
            xiData_radians = np.radians(xiData_degrees)
            etaData_radians = np.radians(etaData_degrees)
            inputRaster.setGridData(xiData_radians, etaData_radians, variableData)
            # Make sure a surface radius is given
            if surfaceRadius==None:
                msg = 'A surface radius must be specified. When calculating a'+\
                      ' raster gradient magnitude on EPSG:4326 the curvature of'+\
                      ' the surface is taken into account, so you must supply'+\
                      ' the surfaceRadius argument.\n'
                raise Exception(msg)
        elif targetCRS_string != self.coordRefSystem.authid():
            inputRaster.changeCoordRefSystem(targetCRS_string)
        # Calculate raster partial derivatives
        xi_derivative_raster = inputRaster.xiPartialDerivative()
        eta_derivative_raster = inputRaster.etaPartialDerivative()
        gradientMagn = np.zeros_like(variableData)
        for etaCounter in range(numb_etaPoints):
            for xiCounter in range(numb_xiPoints):
                if targetCRS_string == 'EPSG:4326':
                    # Calculate Cartesian Coordinates of current point
                    longitude = xiData_radians[etaCounter,xiCounter]
                    latitude = etaData_radians[etaCounter,xiCounter]
                    x = surfaceRadius*np.cos(latitude)*np.cos(longitude)
                    y = surfaceRadius*np.cos(latitude)*np.sin(longitude)
                    z = surfaceRadius*np.sin(latitude)
                    # Calculate Contra-variant basis vectors, Ei_j = d ksi_i/d x_j
                    # Contra-variant basis vector E1 components:
                    E1_1 = -y/(x**2 + y**2)
                    E1_2 = x/(x**2 + y**2)
                    E1_3 = 0.0
                    # Contra-variant basis vector E2 components:
                    E2_1 = -x*z/((x**2 + y**2 + z**2)*np.sqrt(x**2 + y**2))
                    E2_2 = -y*z/((x**2 + y**2 + z**2)*np.sqrt(x**2 + y**2))
                    E2_3 = (x**2 + y**2)/((x**2 + y**2 + z**2)*np.sqrt(x**2 + y**2))
                    # Calculate gradient vector at current point, using
                    # the contra-variant basis vectors.
                    gradientVectorComponents = [
                      xi_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E1_1 + \
                      eta_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E2_1,
                      xi_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E1_2 + \
                      eta_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E2_2,
                      xi_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E1_3 + \
                      eta_derivative_raster.atMatrixIndices(etaCounter, xiCounter)*E2_3
                                               ]
                else:
                    # Calculate gradient vector at current point, using
                    # the scalar partial derivatives only.
                    gradientVectorComponents = [
                      xi_derivative_raster.atMatrixIndices(etaCounter, xiCounter),
                      eta_derivative_raster.atMatrixIndices(etaCounter, xiCounter)
                                               ]
                # Calculate gradient vector magnitude at current point
                gradientMagn[etaCounter,xiCounter] = \
                    np.sqrt(np.sum(np.power(gradientVectorComponents,2.0)))
        resultRaster = raster()
        crs = qgis.core.QgsCoordinateReferenceSystem()
        crs.createFromString(targetCRS_string)
        resultRaster.setCoordRefSystem(crs)
        resultRaster.setGridData(xiData_degrees, etaData_degrees, gradientMagn)
        return resultRaster

    def xiPartialDerivative(self):
        """ Calculate the partial derivative of the variable data with respect to xi. If global, coordinates must be in degrees. """
        #TODO: interrogate CRS to verify coordinates are in degrees when global.
        numb_xiPoints, numb_etaPoints, delta_xi, delta_eta = self.getRasterResolution()
        variable_data = self.getVariableData()
        difference_data = np.gradient(variable_data, delta_xi)
        difference_data_ksi = difference_data[1]
        # The numpy gradient operator uses centre-differences in the interior and
        # one-sided differences in the periphery. If the raster is of global coverage
        # we exploit the East-West periodicity to use centre-differences throughout.
        if self.isGlobal:
            xiVector, etaVector = self.getCoordinatesData()
            xi_east = xiVector[-1]
            xi_west = xiVector[0]
            delta_xi_periodic = 360.0 - abs(xi_east - xi_west)
            # If the grid extends -180 to 180 (or 0 to 360 etc) delta_xi_periodic
            # will be zero. Otherwise there is a "gap" and we must carefully
            # evaluate the derivative.
            if delta_xi_periodic <= 1e-6:
                difference_data_ksi_east = difference_data_ksi[-1,:]
                difference_data_ksi_west = difference_data_ksi[0,:]
                difference_data_ksi[-1,:] = difference_data_ksi[0,:] =\
                      0.5*(difference_data_ksi_east + difference_data_ksi_west)
            else:
                difference_data_ksi_east = difference_data_ksi[-1,:]
                difference_data_ksi_west = difference_data_ksi[0,:]
                variable_data_east = variable_data[-1,:]
                variable_data_west = variable_data[0,:]
                difference_data_ksi_periodic = \
                      (variable_data_west - variable_data_east)/delta_xi_periodic
                difference_data_ksi[-1,:] = \
                      0.5*(difference_data_ksi_east + difference_data_ksi_periodic)
                difference_data_ksi[0,:] =\
                      0.5*(difference_data_ksi_west + difference_data_ksi_periodic)
        # Create new raster object, set data and return.
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        xiData, etaData, variableData = self.getGridData()
        resultRaster.setGridData(xiData, etaData, difference_data_ksi)
        return resultRaster

    def etaPartialDerivative(self):
        """ Calculate the partial derivative of the variable data with respect to eta. """
        numb_xiPoints, numb_etaPoints, delta_xi, delta_eta = self.getRasterResolution()
        variable_data = self.getVariableData()
        difference_data = np.gradient(variable_data, delta_eta)
        difference_data_eta = difference_data[0]
        # Create new raster object, set data and return.
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        xiData, etaData, variableData = self.getGridData()
        resultRaster.setGridData(xiData, etaData, difference_data_eta)
        return resultRaster

    def atGridIndices(self, xiIndex, etaIndex):
        """ Return the field (raster) value at a given grid point.

        This method returns the raster variable at a given grid point,
        where the grid indices are specified in the following order:
        xi index first, follwed by the eta index. For example longitude
        index, latitude index. Note that the internal storge of
        the variable is different: eta, ksi; see method atMatrixIndices()

        :arg int xiIndex: The xi-index of the grid point.
        :arg int etaIndex: The eta-index of the grid point.
        :returns: The value of the raster variable at the specified point.
        """
        return self.variableData[etaIndex, xiIndex]

    def atMatrixIndices(self, etaIndex, xiIndex):
        """ Return the field (raster) value at a given matrix point.

        This method returns the raster variable at a given matrix point,
        where the indices are specified in the following order:
        eta index first, followed by the xi index. For example, latitude
        index followed by longitude index.

        :arg int etaIndex: The eta-index of the grid point.
        :arg int xiIndex: The xi-index of the grid point.
        :returns: The value of the raster variable at the specified point.
        """
        return self.variableData[etaIndex, xiIndex]

    def __add__(self, other):
        """ Add another raster object to this one.

        This method adds structured field data, point by point. Typically,
        the grids the two fields are defined on must be the same. Here
        different grids can be given as input, but this method will
        re-grid (via interpolation) the second field (right-hand-side
        operand) to match the first (left-hand-side operand).
        The same applies to the coordinate reference systems; If they
        differ, the second is re-projected onto the CRS of the first field.

        :arg qmesh.raster self: Left-hand-side operand
        :arg qmesh.raster other: Right-hand-side operand
        :rtype: qmesh.rasterTools.raster
        :returns: A qmesh.rasterTools.raster object containing the point-by-point addition result. This raster is on the coordinate reference system and grid of the left-hand-side operand.
        """
        # Check both rasters have same CRS. If not, change the CRS of the
        # second raster, to match the CRS of the first.
        self_coordRefSystem = self.getCoordRefSystem()
        other_coordRefSystem = other.getCoordRefSystem()
        if self_coordRefSystem != other_coordRefSystem :
            other.changeCoordRefSystem(self_coordRefSystem.authid())
        # Check both rasters have same grid. If not, re-grid second
        # field to match the first argument.
        self_xiData, self_etaData, self_variableData = self.getGridData()
        other_xiData, other_etaData, other_variableData = other.getGridData()
        if self_xiData != other_xiData and self_etaData != other_etaData:
            self_xiMin, self_xiMax, self_etaMin, self_etaMax = self.getRasterBounds()
            self_numb_xiPoints, self_numb_etaPoints, self_delta_xi, self_delta_eta = \
                    self.getRasterResolution()
            other.reGrid(self_xiMin, self_xiMax, self_etaMin, self_etaMax, \
                         self_numb_xiPoints, self_numb_etaPoints)
        # Add grid data and return new object
        sum_variableData = self_variableData + other_variableData
        sumRaster = raster()
        sumRaster.setCoordRefSystem(self_coordRefSystem)
        sumRaster.setGridData(self_xiData, self_etaData, sum_variableData)
        return sumRaster

    def constantDataMultiply(self, constant):
        """ Multiply the variable data in the raster object by a constant value.
        
        :arg float constant: The constant value that the data should be multiplied with.
        :return: The resulting raster object containing the new data.
        """
        self_xiData, self_etaData, self_variableData = self.getGridData()
        result_variableData = self_variableData * constant
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        resultRaster.setGridData(self_xiData, self_etaData, result_variableData)
        return resultRaster

    def constantDataAdd(self, constant):
        """ Add a constant value onto the variable data in the raster object.
        
        :arg float constant: The constant value to be added.
        :return: The resulting raster object containing the new data.
        """
        self_xiData, self_etaData, self_variableData = self.getGridData()
        result_variableData = self_variableData + constant
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        resultRaster.setGridData(self_xiData, self_etaData, result_variableData)
        return resultRaster

    def linearDataTransform(self, outputValue1=None, outputValue2=None,
                            inputValue1=None, inputValue2=None,
                            transformAlpha=None, transformBeta=None):
        """ Apply a linear transformation to the raster data.

        This method calculates the linear transformation alpha*z + beta
        on the raster object is invoked on and returns the result as a new
        raster object. z is the data on the grid of the raster, the output
        raster is defined on the same grid coordinates as the input one. The 
        mapping parameters (alpha & beta) can be calculated from specified
        values of the output and input raster data, alternatively alpha
        and beta can be specified. Note that only one lot (values or
        alpha-&-beta) should be specified. Parameters otherwise specified
        are considered to be erroneous. However, if invoked with no arguments
        the method defaults to outputValue1 = 0, outputValue2 = 1 and
        inputValue1 = min(z), inputValue2 = max(z). If only the output values
        are specified, the method defaults to inputValue1 = min(z),
        inputValue2 = max(z). If only inputValues are specified, outputValue1 = 0,
        outputValue2 = 1

        :arg outputValue1: The minimum of the data of the output raster (optional, defaults to None).
        :arg outputValue2: The maximum of the data of the output raster (optional, defaults to None).
        :arg transformAlpha: alpha (slope) of the linear transformation (optional, defaults to None).
        :arg transformBeta: beta (constant) of the linear transformation (optional, defaults to None).
        :rtype: qmesh.rasterTools.raster
        :returns: A qmesh.rasterTools.raster object containing the linear data transformation.
        :raises Exception: if input arguments are incorrectly specified.
        """
        self_xiData, self_etaData, self_variableData = self.getGridData()
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        # Calculation in case min and max are specified
        if outputValue1!=None and outputValue2!=None and \
           inputValue1==None and inputValue2==None and \
           transformAlpha==None and transformBeta==None:
            inputValue1, inputValue2 = self.getDataBounds()
            transformAlpha = (outputValue2 - outputValue1)/(inputValue2 - inputValue1)
            transformBeta = outputValue1 - transformAlpha*inputValue1
            result_variableData = transformAlpha*self_variableData + transformBeta
            resultRaster.setGridData(self.xiData, self.etaData, result_variableData)
        # Calculation in case min and max are specified
        elif outputValue1!=None and outputValue2!=None and \
           inputValue1!=None and inputValue2!=None and \
           transformAlpha==None and transformBeta==None:
            transformAlpha = (outputValue2 - outputValue1)/(inputValue2 - inputValue1)
            transformBeta = outputValue1 - transformAlpha*inputValue1
            result_variableData = transformAlpha*self_variableData + transformBeta
            resultRaster.setGridData(self.xiData, self.etaData, result_variableData)
        # Calculation in case transformation parameters are specified
        elif outputValue1==None and outputValue2==None and \
           inputValue1==None and inputValue2==None and \
           transformAlpha!=None and transformBeta!=None:
            result_variableData = transformAlpha * self_variableData + transformBeta
            resultRaster.setGridData(self.xiData, self.etaData, result_variableData)
        # Default to outputValue1=0.0, outputValue2=1.0
        elif outputValue1==None and outputValue2==None and \
           inputValue1!=None and inputValue2!=None and \
           transformAlpha==None and transformBeta==None:
            outputValue1=0.0
            outputValue2=1.0
            transformAlpha = (outputValue2 - outputValue1)/(inputValue2 - inputValue1)
            transformBeta = outputValue1 - transformAlpha*inputValue1
            result_variableData = transformAlpha*self_variableData + transformBeta
            resultRaster.setGridData(self.xiData, self.etaData, result_variableData)
        # Default to outputValue1=0.0, outputValue2=1.0
        elif outputValue1==None and outputValue2==None and \
           inputValue1==None and inputValue2==None and \
           transformAlpha==None and transformBeta==None:
            outputValue1=0.0
            outputValue2=1.0
            inputValue1, inputValue2 = self.getDataBounds()
            transformAlpha = (outputValue2 - outputValue1)/(inputValue2 - inputValue1)
            transformBeta = outputValue1 - transformAlpha*inputValue1
            result_variableData = transformAlpha*self_variableData + transformBeta
            resultRaster.setGridData(self.xiData, self.etaData, result_variableData)
        # Arguments are incorrectly specified at input
        else:
            msg = 'ERROR: Input arguments are incorrect. You must'+\
                  ' either specify min & max or alpha & beta of the'+\
                  ' transformation. The arguments are:\n'+\
                  '   outputValue1: '+str(outputValue1)+'\n'+\
                  '   outputValue2: '+str(outputValue2)+'\n'+\
                  '   transformAlpha: '+str(transformAlpha)+'\n'+\
                  '   transformBeta: '+str(transformBeta)+'\n'
            raise Exception(msg)
        return resultRaster

    def __abs__(self):
        """ Take the absolute value of the variable data.
        
        :returns: The resulting raster object.
        """
        self_xiData, self_etaData, self_variableData = self.getGridData()
        result_variableData = np.abs(self_variableData)
        resultRaster = raster()
        resultRaster.setCoordRefSystem(self.getCoordRefSystem())
        resultRaster.setGridData(self_xiData, self_etaData, result_variableData)
        return resultRaster
        
    def boundMinimum(self, minimumValue):
        """ Bound the values of the variable data below by some specified value.
        
        :arg float minimumValue: The minimum value that the variable data is allowed to take.
        """
        etaIindices, xiIndices = np.where(self.getVariableData() <= minimumValue)
        for etaCounter, xiCounter in zip(etaIindices, xiIndices):
            self.variableData[etaCounter,xiCounter] = minimumValue
            
    def boundMaximum(self, maximumValue):
        """ Bound the values of the variable data above by some specified value.
        
        :arg float maximumValue: The maximum value that the variable data is allowed to take.
        """
        etaIindices, xiIndices = np.where(self.getVariableData() >= maximumValue)
        for etaCounter, xiCounter in zip(etaIindices, xiIndices):
            self.variableData[etaCounter,xiCounter] = maximumValue
        

def maximumRaster(rasterList, ignore_nans=True):
    """ Calculate the point-by-point maximum between multiple
    rasters. Typically, the grids the fields are defined on must be the
    same. Here different grids can be given as input, but this method will
    re-grid (via interpolation) the fields to match the first.
    The same applies to the coordinate reference systems; If they
    differ, they are re-projected onto the CRS of the first field.

    :arg qmesh.raster rasterList: A list of rasters. Must contain more than one raster.
    :rtype: qmesh.raster
    :returns: A qmesh.raster object containing the point-by-point minimum across all rasters. This raster is on the coordinate reference system and grid jof the left-hand-side operand.
    :raises Exception: if the input list holds less than 2 objects.
    """
    return minmax(rasterList, operator='max', ignore_nans=ignore_nans)

def minimumRaster(rasterList, ignore_nans=True):
    """ Calculate the point-by-point minimum between multiple
    rasters. Typically, the grids the fields are defined on must be the
    same. Here different grids can be given as input, but this method will
    re-grid (via interpolation) the fields to match the first.
    The same applies to the coordinate reference systems; If they
    differ, they are re-projected onto the CRS of the first field.

    :arg qmesh.raster rasterList: A list of rasters. Must contain more than one raster.
    :rtype: qmesh.raster
    :returns: A qmesh.raster object containing the point-by-point minimum across all rasters. This raster is on the coordinate reference system and grid jof the left-hand-side operand.
    :raises Exception: if the input list holds less than 2 objects.
    """
    return minmax(rasterList, operator='min', ignore_nans=ignore_nans)

def minmax(rasterList, operator=None, ignore_nans=True):
    """ Calculate the point-by-point minimum or maximum between multiple
    rasters. Typically, the grids the fields are defined on must be the
    same. Here different grids can be given as input, but this method will
    re-grid (via interpolation) the fields to match the first.
    The same applies to the coordinate reference systems; If they
    differ, they are re-projected onto the CRS of the first field.

    :arg qmesh.raster rasterList: A list of rasters. Must contain more than one raster.
    :rtype: qmesh.raster
    :returns: A qmesh.raster object containing the point-by-point minimum across all rasters. This raster is on the coordinate reference system and grid jof the left-hand-side operand.
    :raises Exception: if the input list holds less than 2 objects.
    """
    # Check input list contains more than 1
    if len(rasterList) < 2:
        msg = 'ERROR: List of rasters must contain more than'+\
              ' one raster.\n'
        raise Exception(msg)
    # Check all rasters have same CRS
    firstRaster = rasterList.pop(0)
    coordRefSystem = firstRaster.getCoordRefSystem()
    for nextRaster in rasterList:
        other_coordRefSystem = nextRaster.getCoordRefSystem()
        if coordRefSystem != other_coordRefSystem :
            nextRaster.changeCoordRefSystem(coordRefSystem.authid())
    # Check rasters have same grid. If not, re-grid
    # field to match the first one.
    xiData, etaData, variableData = firstRaster.getGridData()
    for nextRaster in rasterList:
        other_xiData, other_etaData, other_variableData = nextRaster.getGridData()
        if np.any(xiData != other_xiData) and np.any(etaData != other_etaData):
            xiMin, xiMax, etaMin, etaMax = firstRaster.getRasterBounds()
            numb_xiPoints, numb_etaPoints, delta_xi, delta_eta = \
                    firstRaster.getRasterResolution()
            nextRaster.reGrid(xiMin, xiMax, etaMin, etaMax, \
                         numb_xiPoints, numb_etaPoints)
            other_xiData, other_etaData, other_variableData = nextRaster.getGridData()
        # Accumulate min or max
        if operator=='min' and not ignore_nans:
            variableData = np.minimum(variableData, other_variableData)
        elif operator=='min' and ignore_nans:
            variableData = np.fmin(variableData, other_variableData)
        elif operator=='max' and not ignore_nans:
            variableData = np.maximum(variableData, other_variableData)
        elif operator=='max' and ignore_nans:
            variableData = np.fmax(variableData, other_variableData)
        else:
            raise Exception('operator must be min or max')
    # Set-up and return new object
    minmaxRaster = raster()
    minmaxRaster.setCoordRefSystem(coordRefSystem)
    minmaxRaster.setGridData(xiData, etaData, variableData)
    return minmaxRaster


class gradationToShapes(raster):
    """ Class for calculating gradations to shapes - lines or polygons. """
    from ..vector.shapefileTools import Shapes
    def __init__(self, shapes = Shapes()):
        """ Class initialisation. Set all the attributes to None, except for the shape objects. """
        super(gradationToShapes, self).__init__()
        self.shapes = shapes
        self.metricAtShapes = None
        self.metricAwayFromShapes = None
        self.gradationDistance = None
        self.coordsys = None
        
    def setShapes(self, shapes):
        """ Set the shapes in the raster object. """
        self.shapes = shapes
        self.coordsys = shapes.getCoordRefSystem()
        
    def setRasterBounds(self, xiMin, xiMax,
                              etaMin, etaMax):
        """ Set the bounds of the raster object.
        
        :arg float xiMin: The minimum value of xi.
        :arg float xiMax: The maximum value of xi.
        :arg float etaMin: The minimum value of eta.
        :arg float etaMax: The maximum value of eta.
        """
        self.xiMin = xiMin
        self.xiMax = xiMax
        self.etaMin = etaMin
        self.etaMax = etaMax
        
    def setRasterResolution(self, numb_xiPoints,
                                  numb_etaPoints):
        """ Set the resolution of the raster.
        
        :arg int numb_xiPoints: The number of xi points.
        :arg int numb_etaPoints: The number of eta points.
        """
        self.numb_xiPoints = numb_xiPoints
        self.numb_etaPoints = numb_etaPoints
        
    def setGradationParameters(self, metricAtShapes,
                               metricAwayFromShapes,
                               gradationDistance,
                               gradationStartDistance=0.0):
        """ Set the gradation parameters. """
        self.metricAtShapes = metricAtShapes 
        self.metricAwayFromShapes = metricAwayFromShapes
        self.gradationDistance = gradationDistance
        self.gradationStartDistance = gradationStartDistance

    def writeDistanceNetCDF(self, outputRasterFilename):
        """ Write the distance function in NetCDF format.
        
        :arg str outputRasterFilename: The desired name of the output file.
        """
        #TODO: try other netCDF interfaces.
        from datetime import datetime
        import os
        import subprocess
        from ..config import _subprocess_log_queue
        # Construct temporary file-name, for storing shapes as well as "rasterised" shapes
        time = datetime.now()
        shapesFilenamePattern = '/tmp/shapes'+time.isoformat()
        shapesSHPFilename = '/tmp/shapes'+time.isoformat()+'.shp'
        rasterisedShapesFilename = '/tmp/rasterisedShapes'+time.isoformat()+'.nc'
        # Calculate output resolution
        delta_xi = abs(float(self.xiMax) - float(self.xiMin))/(float(self.numb_xiPoints)-1.)
        delta_eta = abs(float(self.etaMax) - float(self.etaMin))/(float(self.numb_etaPoints)-1.)
        # The distance function calculation is here carried out using gdal_proximity.
        # In order for that utility to give accurate distances the pixels must be
        # square: delta_xi == delta_eta. Below we check for that condition, 
        # and if not the case we change the resolution and/or enlarge the raster
        # extents, to make the pixels square.
        if delta_xi != delta_eta:
            LOG.info('The distance function raster grid must have square cells (GDAL requirement).')
            LOG.info(' This is not possible with the specified raster-resolution.')
            if delta_xi < delta_eta:
                delta_eta = delta_xi
                numb_etaPoints = (abs(self.etaMax - self.etaMin)/delta_eta) + 1
                self.numb_etaPoints = int(np.ceil(numb_etaPoints))
                self.etaMax = delta_eta*(self.numb_etaPoints - 1) + self.etaMin
            elif delta_eta < delta_xi:
                delta_xi = delta_eta
                numb_xiPoints = (abs(self.xiMax - self.xiMin)/delta_xi) + 1
                self.numb_xiPoints = int(np.ceil(numb_xiPoints))
                self.xiMax = delta_xi*(self.numb_xiPoints - 1) + self.xiMin
            LOG.info(' Changed raster extents and resolution to:')
            LOG.info(' '+str(self.numb_xiPoints)+' points along first direction, extending from '+str(self.xiMin)+' to '+str(self.xiMax)+'.')
            LOG.info(' '+str(self.numb_etaPoints)+' points along second direction, extending from '+str(self.etaMin)+' to '+str(self.etaMax)+'.')
        # Delete temporary files if they already exist
        filesForDeletion = [shapesFilenamePattern+'.shp',
                            shapesFilenamePattern+'.shx',
                            shapesFilenamePattern+'.qpj',
                            shapesFilenamePattern+'.prj',
                            shapesFilenamePattern+'.cpg',
                            shapesFilenamePattern+'.dbf',
                            rasterisedShapesFilename]
        for fileName in filesForDeletion:
            try:
                os.remove(fileName)
            except OSError: #in case given file does not exist
                pass
        # Write shapes file, but do not output to log file
        logLevel = LOG.level
        LOG.setLevel('WARNING')
        self.shapes.writeFile(shapesSHPFilename)
        LOG.setLevel(logLevel)
        # "Rasterize" the shapes
        qmesh_pid = subprocess.os.getpid()
        gdal_rasterize_stdoutFileName = '/tmp/gdal_rasterize_stdout_'+str(qmesh_pid)
        gdal_rasterize_stdout = open(gdal_rasterize_stdoutFileName,'w')
        gdal_rasterize_stderrFileName = '/tmp/gdal_rasterize_stderr_'+str(qmesh_pid)
        gdal_rasterize_stderr = open(gdal_rasterize_stderrFileName,'w')
        gdal_rasterize_cmd = ['gdal_rasterize',
              '-q','-burn','1','-a_nodata','0','-init','0','-at',
              '-tr',str(delta_xi),str(delta_eta),
              '-te',str(self.xiMin),str(self.etaMin),str(self.xiMax),str(self.etaMax),
              '-of','netCDF',shapesSHPFilename,rasterisedShapesFilename]
        gdal_rasterize_proc = \
          subprocess.Popen(gdal_rasterize_cmd,
                           stdout=gdal_rasterize_stdout,
                           stderr=gdal_rasterize_stderr)
        gdal_rasterize_log_queue = _subprocess_log_queue(gdal_rasterize_proc,
                                                         'gdal_rasterize',
                                                         gdal_rasterize_stdoutFileName,
                                                         gdal_rasterize_stderrFileName)
        #Get a child logger (child of the qmesh logger) to output stdout and stderr from gdal_proximity
        gdal_rasterize_logger = logging.getLogger('qmesh.mesh.gdal_rasterize')
        while True:
          if gdal_rasterize_log_queue.emptyQueue() and gdal_rasterize_proc.poll() is not None:
            break
          elif gdal_rasterize_log_queue.emptyQueue() and gdal_rasterize_proc.poll() is None:
            continue
          elif not gdal_rasterize_log_queue.emptyQueue():
            line = gdal_rasterize_log_queue.getLine()
            log_message = line.strip('\n').split()
            if log_message[0] == 'Warning':
              gdal_rasterize_logger.warning(' '.join(log_message[2:]))
            elif log_message[0] == 'ERROR':
              gdal_rasterize_logger.error(' '.join(log_message[2:]))
            else:
              gdal_rasterize_logger.info(' '.join(log_message))
        try:
            if gdal_rasterize_proc.returncode != 0:
              msg='GDAL rasterize error'
              raise Exception(msg)
        except:
            LOG.error(msg, exc_info=True)
            raise
        # Delete distance-function raster file if it already exists
        try:
            os.remove(outputRasterFilename)
        except OSError: #in case given file does not exist
            pass
        # close file handles
        gdal_rasterize_stderr.close()
        gdal_rasterize_stdout.close()

        # Distance-function calulation, outputing NetCDF file.
        gdal_proximity_cmd = ['gdal_proximity.py',
                              '-q',
                              rasterisedShapesFilename,
                              outputRasterFilename,
                              '-of', 'netCDF', '-distunits', 'GEO']
        gdal_proximity_stdoutFileName = '/tmp/gdal_proximity_stdout_'+str(qmesh_pid)
        gdal_proximity_stdout = open(gdal_proximity_stdoutFileName,'w')
        gdal_proximity_stderrFileName = '/tmp/gdal_proximity_stderr_'+str(qmesh_pid)
        gdal_proximity_stderr = open(gdal_proximity_stderrFileName,'w')
        gdal_proximity_proc = \
          subprocess.Popen(gdal_proximity_cmd,
                           stdout=gdal_proximity_stdout,
                           stderr=gdal_proximity_stderr)
        gdal_proximity_log_queue = _subprocess_log_queue(gdal_proximity_proc,
                                                         'gdal_proximity',
                                                         gdal_proximity_stdoutFileName,
                                                         gdal_proximity_stderrFileName)
        #Get a child logger (child of the qmesh logger) to output stdout and stderr from gdal_proximity
        gdal_proximity_logger = logging.getLogger('qmesh.mesh.gdal_proximity')
        while True:
          if gdal_proximity_log_queue.emptyQueue() and gdal_proximity_proc.poll() is not None:
            break
          elif gdal_proximity_log_queue.emptyQueue() and gdal_proximity_proc.poll() is None:
            continue
          elif not gdal_proximity_log_queue.emptyQueue():
            line = gdal_proximity_log_queue.getLine()
            log_message = line.strip('\n').split()
            if log_message[0] == 'Warning':
              gdal_proximity_logger.warning(' '.join(log_message[2:]))
            elif log_message[0] == 'ERROR':
              gdal_proximity_logger.error(' '.join(log_message[2:]))
            else:
              gdal_proximity_logger.info(' '.join(log_message))
        try:
            if gdal_proximity_proc.returncode !=0 :
              msg='GDAL proximity error'
              raise Exception()
        except:
            LOG.error(msg, exc_info=True)
            raise
        # clea up file handles
        gdal_proximity_stderr.close()
        gdal_proximity_stdout.close()


        # Clean-up temporary files.
        filesForDeletion = [shapesFilenamePattern+'.shp',
                            shapesFilenamePattern+'.shx',
                            shapesFilenamePattern+'.qpj',
                            shapesFilenamePattern+'.prj',
                            shapesFilenamePattern+'.cpg',
                            shapesFilenamePattern+'.dbf',
                            rasterisedShapesFilename]
        for fileName in filesForDeletion:
            try:
                os.remove(fileName)
            except OSError: #in case given file does not exist
                pass


    def calculateLinearGradation(self):
        """ Compute the linear gradation. """
        from datetime import datetime
        import os
        LOG.info('Calculating linear gradation from shapes')
        LOG.debug('    Gradating from '+str(self.metricAtShapes)+' to '+str(self.metricAwayFromShapes))
        # Construct temporary file-names.
        time = datetime.now()
        rasterisedDistanceFile = '/tmp/rasterisedDistance'+time.isoformat()+'.nc'
        # Calculate and read-in distance-function.
        self.writeDistanceNetCDF(rasterisedDistanceFile)
        logLevel = LOG.level
        LOG.setLevel('WARNING')
        # this can reset the coord ref system to None or ""
        self.fromFile(rasterisedDistanceFile)
        LOG.setLevel(logLevel)
        # Map distance to mesh-size-metric
        meshSizeMetric = self.variableData*((float(self.metricAwayFromShapes) - float(self.metricAtShapes))/float(self.gradationDistance)) + \
              float(self.metricAwayFromShapes) - (float(self.gradationDistance) + float(self.gradationStartDistance))*\
              ((float(self.metricAwayFromShapes) - float(self.metricAtShapes))/float(self.gradationDistance))
        meshSizeMetric = np.fmax(meshSizeMetric, self.metricAtShapes)
        meshSizeMetric = np.fmin(meshSizeMetric, self.metricAwayFromShapes)
        # force reset of coord system
        self.coordRefSystem = self.coordsys
        self.variableData = meshSizeMetric
        # Clean-up temporary files.
        filesForDeletion = [rasterisedDistanceFile, rasterisedDistanceFile+'.aux.xml']
        for fileName in filesForDeletion:
            try:
                os.remove(fileName)
            except OSError: # In case given file does not exist
                pass

def elementSizeWaveCFL(bathymetryRaster, targetRatio=100.0, wavePeriod=44712.0, gravityAcceleration=9.81):
    """ Calculate the element size based on the Courant number. """
    resultRaster = bathymetryRaster.deepcopy() 
    resultData = (np.sqrt(np.abs(resultRaster.getVariableData())*gravityAcceleration)/targetRatio)*wavePeriod
    resultRaster.setVariableData(resultData)
    return resultRaster
