#!/usr/bin/env python


def OrkneyShetlandIsles_UTM30_InnerSoundTurbines():
    '''Todo: add docstring '''
    import numpy as np
    import qgis.core
    import sys
    import qmesh
    #Initialising qgis API
    qmesh.initialise()
    #Reading in the shapefile describing the domain boundaries, and creating a gmsh file.
    boundaries = qmesh.vector.Shapes()
    boundaries.fromFile('OrkneyShetlandIsles_singleRegion.shp')
    loopShapes = qmesh.vector.identifyLoops(boundaries,
              isGlobal=False, defaultPhysID=1000,
              fixOpenLoops=True)
    polygonShapes = qmesh.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=50000, 
                                                                     meshedAreaPhysID = 1)
    #Create loop and polygon for Inner Sound.
    innerSound_plot_lines = qmesh.vector.Shapes()
    innerSound_plot_lines.fromFile('innerSound_plot_lines.shp')
    innerSound_plot_loops = qmesh.vector.identifyLoops(innerSound_plot_lines,
              fixOpenLoops=True, extraPointsPerVertex=10)
    innerSound_plot_polygon = qmesh.vector.identifyPolygons(innerSound_plot_loops, 
                                                                     meshedAreaPhysID = 2)
    #Using a random-number generator to distribute points inside
    # the Inner Sound plot, with the additional constrains that
    # they must be 200m apart. Locate 100 such points. Then create
    # a cicle for each point, representative of a tidal turbine.
    #Seed the random-number generator so as to obtain same set of
    # points every time the test is run. 
    np.random.seed(1)
    #Calculate mapping
    # parameters so as to map the interval [0,1] to an area
    # enclosing the tidal plot. The numbers used to form
    # the parameters below are corner
    # coordinates from a rectangle -in UTM30- that ecloses
    # the tidal plot
    ksi_param_1 = 494800 - 490400
    ksi_param_2 = 490400
    eta_param_1 = 6503800 - 6501600
    eta_param_2 = 6501600
    turbinePoints = []
    #Change the coordinate rreference system of the Inner Sound plot
    # from lon-lat into UTM30 .
    innerSound_plot_polygon.changeCoordRefSystem('EPSG:32630')
    polygon = innerSound_plot_polygon.getFeatures()[0]
    sys.stdout.write('Placing the turbines...\n')
    sys.stdout.flush()
    while len(turbinePoints) < 180:
        #Get new "arbitary" point, carry out mapping using
        # the parameters calulated above.
        newPointX = np.random.random()*ksi_param_1 + ksi_param_2
        newPointY = np.random.random()*eta_param_1 + eta_param_2
        #Construct a large circle around that point, 20m in radius
        # Then check if it is included into the tidal plot. The larger
        # radious should eliminate the possibility of the centre-point
        # of the cicle lying inside the polygon, but too close to the
        # plot boundary such that a part of the tidal-turbine is
        # outside the plot.
        turbineBigCircle = qmesh.vector.primitiveShapes.Circle(newPointX, newPointY, 15.0, 20, 'EPSG:32630')
        turbineBigPolygon = turbineBigCircle.asQgsPolygonFeature()
        if polygon.geometry().contains(turbineBigPolygon.geometry()):
            if len(turbinePoints) == 0:
                turbinePoints.append((newPointX, newPointY))
            else:
                isNewPointAdmissible = True
                for otherPoint in turbinePoints:
                    distance = np.sqrt( (otherPoint[0] - newPointX)**2 + (otherPoint[1] - newPointY)**2)
                    if distance < 100.0:
                        isNewPointAdmissible = False
                        break
                if isNewPointAdmissible:
                    turbinePoints.append((newPointX, newPointY))
    #Write the turbine "center" points into an asci file, usefull for declaring detectors
    # in the flml
    np.savetxt('turbineDetectors_UTM30N.txt', turbinePoints, header='#x,y', comments='#UTM30N\n')
    #Write the turbines as circular areas of 20m diameter.
    turbines = qmesh.vector.Shapes()
    turbines.setCoordRefSystemFromString('EPSG:32630')
    turbines.setShapeType(qgis.core.QgsWkbTypes.LineString)
    for turbinePoint in turbinePoints:
        turbineCircle = qmesh.vector.primitiveShapes.Circle(turbinePoint[0], turbinePoint[1], 10.0, 20, 'EPSG:32630')
        turbines.addFeature(turbineCircle.asQgsFeature())
    turbines.writeFile('turbines.shp')
    turbines.changeCoordRefSystem('EPSG:4326')
    turbineLoops = qmesh.vector.identifyLoops(turbines,
              fixOpenLoops=True)
    turbinePolygons = qmesh.vector.identifyPolygons(turbineLoops,
              meshedAreaPhysID = 3)
    #Reconstruct innerSound_plot_polygon object, so as to obtain
    # it back to EPSG:4326, recall we changed its CRS to EPSG:32630
    # While we could invoke the changeCoordRefSystem() method on it
    # we will not obtain the same point coordinates as before, leading
    # to problems later on.
    innerSound_plot_polygon = qmesh.vector.identifyPolygons(innerSound_plot_loops, 
                                                                     meshedAreaPhysID = 2)
    #Insert tidal plots and turbines into domain.
    domainLines, domainPolygons = qmesh.vector.insertRegions(
                                          loopShapes, polygonShapes,\
                                          innerSound_plot_loops, innerSound_plot_polygon)
    domainLines, domainPolygons = qmesh.vector.insertRegions(
                                          domainLines, domainPolygons,\
                                          turbineLoops, turbinePolygons)
    domainPolygons.writeFile('domainPolygons')
    #Create raster for mesh gradation towards full-resolution shorelines.
    GSHHS_fine_boundaries = qmesh.vector.Shapes()
    GSHHS_fine_boundaries.fromFile('GSHHS_f_L1_lines.shp')
    gradationRaster_shoreline = qmesh.raster.meshMetricTools.gradationToShapes()
    gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
    gradationRaster_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
    gradationRaster_shoreline.setRasterResolution(1500,1500)
    gradationRaster_shoreline.setGradationParameters(150.0,15000.0,1.0)
    gradationRaster_shoreline.calculateLinearGradation()
    gradationRaster_shoreline.writeNetCDF('gradation_to_GSHHS_f_lines.nc')
    #Create raster for mesh gradation towards Shetlands shorelines.
    shetlands_shorelines = qmesh.vector.Shapes()
    shetlands_shorelines.fromFile('Shetlands_shoreline_gebco08.shp')
    gradationRaster_shetlands_shoreline = qmesh.raster.meshMetricTools.gradationToShapes()
    gradationRaster_shetlands_shoreline.setShapes(shetlands_shorelines)
    gradationRaster_shetlands_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
    gradationRaster_shetlands_shoreline.setRasterResolution(1500,1500)
    gradationRaster_shetlands_shoreline.setGradationParameters(1500.0,15000.0,0.5)
    gradationRaster_shetlands_shoreline.calculateLinearGradation()
    gradationRaster_shetlands_shoreline.writeNetCDF('gradation_to_Shetlands_shoreline.nc')
    #Create raster for mesh gradation towards 0m gebco contour on the Scottish mainland coast
    GEBCO08_0mContour = qmesh.vector.Shapes()
    GEBCO08_0mContour.fromFile('GEBCO08_0mContour.shp')
    gradationRaster_GEBCO08_0mContour = qmesh.raster.meshMetricTools.gradationToShapes()
    gradationRaster_GEBCO08_0mContour.setShapes(GEBCO08_0mContour)
    gradationRaster_GEBCO08_0mContour.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
    gradationRaster_GEBCO08_0mContour.setRasterResolution(1500,1500)
    gradationRaster_GEBCO08_0mContour.setGradationParameters(1500.0,15000.0,0.5)
    gradationRaster_GEBCO08_0mContour.calculateLinearGradation()
    gradationRaster_GEBCO08_0mContour.writeNetCDF('gradationRaster_GEBCO08_0mContour.nc')
    #Create raster for mesh gradation towards GSHHS high-resolution lines on the Scottish mainland coast
    GSHHS_h_L1_lines = qmesh.vector.Shapes()
    GSHHS_h_L1_lines.fromFile('GSHHS_h_L1_lines.shp')
    gradationRaster_GSHHS_h_L1_lines = qmesh.raster.meshMetricTools.gradationToShapes()
    gradationRaster_GSHHS_h_L1_lines.setShapes(GSHHS_h_L1_lines)
    gradationRaster_GSHHS_h_L1_lines.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
    gradationRaster_GSHHS_h_L1_lines.setRasterResolution(1500,1500)
    gradationRaster_GSHHS_h_L1_lines.setGradationParameters(1500.0,15000.0,0.5)
    gradationRaster_GSHHS_h_L1_lines.calculateLinearGradation()
    gradationRaster_GSHHS_h_L1_lines.writeNetCDF('GSHHS_h_L1_lines_gradationRaster_.nc')
    #Create raster for mesh gradation towards Inner Sound
    gradationRaster_innerSound_plot = qmesh.raster.meshMetricTools.gradationToShapes()
    gradationRaster_innerSound_plot.setShapes(innerSound_plot_polygon)
    gradationRaster_innerSound_plot.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
    gradationRaster_innerSound_plot.setRasterResolution(1500,1500)
    gradationRaster_innerSound_plot.setGradationParameters(5.0,15000.0,1.0,0.01)
    gradationRaster_innerSound_plot.calculateLinearGradation()
    gradationRaster_innerSound_plot.writeNetCDF('gradation_to_InnerSoundPlot.nc')
    #Calculate overall mesh-metric raster
    meshMetricRaster = qmesh.raster.meshMetricTools.minimumRaster([gradationRaster_shoreline, \
                                                                   gradationRaster_shetlands_shoreline, \
                                                                   gradationRaster_GEBCO08_0mContour, \
                                                                   gradationRaster_GSHHS_h_L1_lines, \
                                                                   gradationRaster_innerSound_plot])
    meshMetricRaster.writeNetCDF('meshMetric.nc')
    #Create domain object and write gmsh files.
    domain = qmesh.mesh.Domain()
    domain.setGeometry(domainLines, domainPolygons)
    #domain.setGeometry(loopShapes, polygonShapes)
    domain.setMeshMetricField(meshMetricRaster)
    domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)
    #Meshing
    domain.gmsh(geoFilename='OrkneyShetlandIsles_UTM30_InnerSoundTurbines.geo', \
                fldFilename='OrkneyShetlandIsles_UTM30_InnerSoundTurbines.fld', \
                mshFilename='OrkneyShetlandIsles_UTM30_InnerSoundTurbines.msh', \
               )
 
def convertMesh():
    mesh = qmesh.mesh.Mesh()
    mesh.readGmsh('OrkneyShetlandIsles_UTM30_InnerSoundTurbines.msh', 'EPSG:32630')
    mesh.writeShapefile('OrkneyShetlandIsles_UTM30_InnerSoundTurbines')    


if __name__ == '__main__':
    import qmesh
    #Initialising qgis API
    qmesh.initialise()
    OrkneyShetlandIsles_UTM30_InnerSoundTurbines()
    convertMesh()
