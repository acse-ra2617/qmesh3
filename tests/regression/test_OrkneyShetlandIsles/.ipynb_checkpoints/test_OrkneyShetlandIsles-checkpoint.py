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

class TestOrkneyShetlandIsles(unittest.TestCase):
    '''Todo: add test documentation as class docstring '''
    def setUp(self):
        '''Method setting-up test fixture.

        This method is not aimed to be invoked by the user. Various resources needed and
        produced by the test are defined here:
        1) The path of the testing function and the project root path (qmesh).
        2) The filenames, with absolute paths, of files needed by the test-cases.
        3) The filenames, with absolute paths, of files produced by the test-cases.
        '''
        import os
        import sys
        import glob
        #Define list of output files, empty for now
        self.output_files = []
        #Construct absolute path to this test file, without the test file-name
        self.this_path = os.path.dirname(os.path.realpath(__file__))
        #Construct absolute path to the test-type directory
        self.project_test_type_path = os.path.split(self.this_path)[0]
        #Construct absolute path to the tests directory
        self.project_tests_path = os.path.split(self.project_test_type_path)[0]
        #Construct absolute path to the project root
        self.project_root_path = os.path.split(self.project_tests_path)[0]
        #Define filenames used in this test suite.
        #Filename of vector (shapefile) containing all the boundaries
        # for this test-suite. Annotated with PhysID.
        self.all_boundaries_vector_filename = os.path.join(
            self.this_path, 'OrkneyShetlandIsles_singleRegion.shp')
        #Filename of vector (shapefile) containing the boundaries
        # of the Inner Sound Tidal plot.
        self.innerSound_plot_vector_filename = os.path.join(
            self.this_path, 'innerSound_plot_lines.shp')
        #Filename of vector (shapefile) containing Orkney shorelines
        # extracted from the GSHHS full-resolution dataset.
        self.GSHHS_f_vector_filename = os.path.join(
            self.this_path, 'GSHHS_f_L1_lines.shp')
        #Filename of vector (shapefile) containing shorelines
        # extracted from the GSHHS high-resolution dataset.
        self.GSHHS_h_vector_filename = os.path.join(
            self.this_path, 'GSHHS_h_L1_lines.shp')
        #Filename of vector (shapefile) containing a 0m contour
        # extracted from GEBCO08 around the Shetland Islands.
        self.shetlands_shoreline_vector_filename = os.path.join(
            self.this_path, 'Shetlands_shoreline_gebco08.shp')
        #Filename of vector (shapefile) containing a 0m contour
        # extracted from GEBCO08 around the Scottish mainland
        self.scotlald_shoreline_GEBCO08_vector_filename = os.path.join(
            self.this_path, 'GEBCO08_0mContour.shp')
        #Filename of NetCDF raster containing the element edge length distribution, calculated
        # using the GSHHS full-resolution shorelines.
        self.GSHHS_f_gradation_raster_filename = os.path.join(
            self.this_path, 'gradation_to_GSHHS_f_lines.nc')
        self.output_files.append(self.GSHHS_f_gradation_raster_filename)
        #Filename of NetCDF raster containing the element edge length distribution, calculated
        # using the GSHHS high-resolution shorelines.
        self.GSHHS_h_gradation_raster_filename = os.path.join(
            self.this_path, 'gradation_to_GSHHS_h_lines.nc')
        self.output_files.append(self.GSHHS_h_gradation_raster_filename)
        #Filename of NetCDF raster containing the element edge length distribution, calculated
        # using the Shetland shorelines extracted from GEBCO08.
        self.shetlands_shoreline_gradation_raster = os.path.join(
            self.this_path, 'gradation_to_Shetlands_shoreline.nc')
        self.output_files.append(self.shetlands_shoreline_gradation_raster)
        #Filename of NetCDF raster containing the element edge length distribution, calculated
        # using the Scottish mainland shorelines extracted from GEBCO08.
        self.gradation_GEBCO08_0mContour_filename = os.path.join(
            self.this_path,'gradationRaster_GEBCO08_0mContour.nc')
        self.output_files.append(self.gradation_GEBCO08_0mContour_filename)
        #Filename of NetCDF raster containing the element edge length distribution, calculated
        # using the Inner Sound tidal plot boundaries.
        self.gradation_to_InnerSoundPlot_filename = os.path.join(
            self.this_path, 'gradation_to_InnerSoundPlot.nc')
        self.output_files.append(self.gradation_to_InnerSoundPlot_filename)
        #Filename of NetCDF raster containing the gradation metric for the UTM30 problem.
        self.meshMetric_raster_filename = os.path.join(
            self.this_path, 'OrkneyShetlandIsles_meshMetric.nc')
        self.output_files.append(self.meshMetric_raster_filename)
        #Gmsh and shapefile filenames for UTM30 test-case.
        self.UTM30_geo_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_UTM30.geo')
        self.UTM30_fld_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_UTM30.fld')
        self.UTM30_msh_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_UTM30.msh')
        self.UTM30_shp_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_UTM30')
        self.output_files.extend(glob.glob(self.UTM30_shp_filename+'.*'))
        #Gmsh and shapefile filenames for PCC test-case.
        self.PCC_geo_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC.geo')
        self.PCC_fld_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC.fld')
        self.PCC_msh_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC.msh')
        self.PCC_shp_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC')
        self.output_files.extend(glob.glob(self.PCC_shp_filename+'.*'))
        #Gmsh and shapefile filenames for PCC with tidal farm plot test-case.
        self.PCC_inner_sound_plot_geo_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_plot.geo')
        self.PCC_inner_sound_plot_fld_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_plot.fld')
        self.PCC_inner_sound_plot_msh_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_plot.msh')
        self.PCC_inner_sound_plot_shp_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_plot')
        self.output_files.extend(glob.glob(self.PCC_inner_sound_plot_shp_filename+'.*'))
        #Gmsh and shapefile filenames for PCC with tidal farm plot & turbines test-case.
        self.PCC_inner_sound_turbines_geo_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_turbines.geo')
        self.PCC_inner_sound_turbines_fld_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_turbines.fld')
        self.PCC_inner_sound_turbines_msh_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_turbines.msh')
        self.PCC_inner_sound_turbines_shp_filename = os.path.join(
            self.this_path,'OrkneyShetlandIsles_PCC_inner_sound_turbines')
        self.output_files.extend(glob.glob(self.PCC_inner_sound_turbines_shp_filename+'.*'))
        #Filename of domain lines for composite polygon domain test cases.
        self.domain_lines_filename = os.path.join(
            self.this_path,'domain_lines')
        self.output_files.extend(glob.glob(self.domain_lines_filename+'.*'))
        #Filename of domain polygons for composite polygon domain test cases.
        self.domain_polygons_filename = os.path.join(
            self.this_path,'domain_polygons')
        self.output_files.extend(glob.glob(self.domain_polygons_filename+'.*'))
        #Filename of turbine-representing-disk lines.
        self.turbine_lines_filename = os.path.join(
            self.this_path, 'turbines')
        self.output_files.extend(glob.glob(self.turbine_lines_filename+'.*'))

    def remove_file(self, filename):
        '''Remove files, whithout failing for files that do not exist.

        This method uses os.remove to remove files, but captures the exception
        thrown in case the files do not exist. Exceptions due to other reasons
        are still raised.
        '''
        import os, errno
        try:
            os.remove(filename)
        except OSError as error:
            if error.errno != errno.ENOENT:
                raise

    def tearDown(self):
        '''Method taking down the test fixture, tiding up.

        Test-cases typically generate files that can obscure working. This method
        will remove all files identified in the setUp method as output files. Note that
        the output files are ignored by the git set-up.
        '''
        for filename in self.output_files:
            self.remove_file(filename)

    def test_OrkneyShetlandIsles_UTM30(self):
        '''Meshing around the Orkney and Shteland Isles in UTM30 as a regression test.'''
        try:
            #Try importing qmesh. If an ImportError exception is thrown, then
            # the test might be run in a development setting, and see if qmesh
            # can be found at the project root.
            try:
                import qmesh3
            except ImportError:
                sys.path.append(self.project_root_path)
                import qmesh3
            #Initialise qmesh
            qmesh3.initialise()
            #Read-in the shapefile describing the domain boundaries, and creating a gmsh file.
            boundaries = qmesh3.vector.Shapes()
            boundaries.fromFile(self.all_boundaries_vector_filename)
            loopShapes = qmesh3.vector.identifyLoops(boundaries, isGlobal=False,
                                                    defaultPhysID=1000, fixOpenLoops=True)
            polygonShapes = qmesh3.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=50000)
            #Create loops and polygons for Inner Sound
            innerSound_plot_lines = qmesh3.vector.Shapes()
            innerSound_plot_lines.fromFile(self.innerSound_plot_vector_filename)
            innerSound_plot_loops = qmesh3.vector.identifyLoops(innerSound_plot_lines,
                                                               fixOpenLoops=True)
            innerSound_plot_polygon = qmesh3.vector.identifyPolygons(innerSound_plot_loops, 
                                                                    meshedAreaPhysID = 2)
            #Create raster for mesh gradation towards full-resolution shorelines.
            GSHHS_fine_boundaries = qmesh3.vector.Shapes()
            GSHHS_fine_boundaries.fromFile(self.GSHHS_f_vector_filename)
            gradationRaster_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
            gradationRaster_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shoreline.setRasterResolution(800,800)
            gradationRaster_shoreline.setGradationParameters(150.0,15000.0,1.0)
            gradationRaster_shoreline.calculateLinearGradation()
            gradationRaster_shoreline.writeNetCDF(self.GSHHS_f_gradation_raster_filename)
            #Create raster for mesh gradation towards Shetlands shorelines.
            shetlands_shorelines = qmesh3.vector.Shapes()
            shetlands_shorelines.fromFile(self.shetlands_shoreline_vector_filename)
            gradationRaster_shetlands_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shetlands_shoreline.setShapes(shetlands_shorelines)
            gradationRaster_shetlands_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shetlands_shoreline.setRasterResolution(800,800)
            gradationRaster_shetlands_shoreline.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_shetlands_shoreline.calculateLinearGradation()
            gradationRaster_shetlands_shoreline.writeNetCDF(
                self.shetlands_shoreline_gradation_raster)
            #Create raster for mesh gradation towards 0m gebco contour on the Scottish mainland.
            scotlald_shoreline_GEBCO08 = qmesh3.vector.Shapes()
            scotlald_shoreline_GEBCO08.fromFile(self.scotlald_shoreline_GEBCO08_vector_filename)
            gradationRaster_GEBCO08_0mContour = qmesh3.raster.gradationToShapes()
            gradationRaster_GEBCO08_0mContour.setShapes(scotlald_shoreline_GEBCO08)
            gradationRaster_GEBCO08_0mContour.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GEBCO08_0mContour.setRasterResolution(800,800)
            gradationRaster_GEBCO08_0mContour.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GEBCO08_0mContour.calculateLinearGradation()
            gradationRaster_GEBCO08_0mContour.writeNetCDF(self.gradation_GEBCO08_0mContour_filename)
            #Create raster for mesh gradation towards GSHHS high-resolution lines on the Scottish mainland coast
            GSHHS_h_L1_lines = qmesh3.vector.Shapes()
            GSHHS_h_L1_lines.fromFile(self.GSHHS_h_vector_filename)
            gradationRaster_GSHHS_h_L1_lines = qmesh3.raster.gradationToShapes()
            gradationRaster_GSHHS_h_L1_lines.setShapes(GSHHS_h_L1_lines)
            gradationRaster_GSHHS_h_L1_lines.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GSHHS_h_L1_lines.setRasterResolution(800,800)
            gradationRaster_GSHHS_h_L1_lines.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GSHHS_h_L1_lines.calculateLinearGradation()
            gradationRaster_GSHHS_h_L1_lines.writeNetCDF(self.GSHHS_h_gradation_raster_filename)
            #Create raster for mesh gradation towards Inner Sound plot
            gradationRaster_innerSound_plot = qmesh3.raster.gradationToShapes()
            gradationRaster_innerSound_plot.setShapes(innerSound_plot_polygon)
            gradationRaster_innerSound_plot.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_innerSound_plot.setRasterResolution(800,800)
            gradationRaster_innerSound_plot.setGradationParameters(10.0,15000.0,1.0,0.01)
            gradationRaster_innerSound_plot.calculateLinearGradation()
            gradationRaster_innerSound_plot.writeNetCDF(self.gradation_to_InnerSoundPlot_filename)
            #Calculate overall mesh-metric raster
            meshMetricRaster = qmesh3.raster.minimumRaster([gradationRaster_shoreline, \
                                                           gradationRaster_shetlands_shoreline, \
                                                           gradationRaster_GEBCO08_0mContour, \
                                                           gradationRaster_GSHHS_h_L1_lines, \
                                                           gradationRaster_innerSound_plot])
            meshMetricRaster.writeNetCDF(self.meshMetric_raster_filename)
            #Generate mesh
            domain = qmesh3.mesh.Domain()
            domain.setGeometry(loopShapes, polygonShapes)
            domain.setMeshMetricField(meshMetricRaster)
            domain.setTargetCoordRefSystem('EPSG:32630', fldFillValue=1000.0)
            domain.gmsh(geoFilename = self.UTM30_geo_filename, \
                        fldFilename = self.UTM30_fld_filename, \
                        mshFilename = self.UTM30_msh_filename)
        except AssertionError:
            self.assertTrue(False)

    def test_OrkneyShetlandIsles_PCC(self):
        '''Todo: add docstring '''
        try:
            #Try importing qmesh. If an ImportError exception is thrown, then
            # the test might be run in a development setting, and see if qmesh
            # can be found at the project root.
            try:
                import qmesh3
            except ImportError:
                sys.path.append(self.project_root_path)
                import qmesh3
            #Initialise qmesh
            qmesh3.initialise()
            #Read-in the shapefile describing the domain boundaries, and creating a gmsh file.
            boundaries = qmesh3.vector.Shapes()
            boundaries.fromFile(self.all_boundaries_vector_filename)
            loopShapes = qmesh3.vector.identifyLoops(boundaries, isGlobal=False,
                                                    defaultPhysID=1000, fixOpenLoops=True)
            polygonShapes = qmesh3.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=50000)
            #Create loops and polygons for Inner Sound
            innerSound_plot_lines = qmesh3.vector.Shapes()
            innerSound_plot_lines.fromFile(self.innerSound_plot_vector_filename)
            innerSound_plot_loops = qmesh3.vector.identifyLoops(innerSound_plot_lines,
                                                               fixOpenLoops=True)
            innerSound_plot_polygon = qmesh3.vector.identifyPolygons(innerSound_plot_loops, 
                                                                    meshedAreaPhysID = 2)
            #Create raster for mesh gradation towards full-resolution shorelines.
            GSHHS_fine_boundaries = qmesh3.vector.Shapes()
            GSHHS_fine_boundaries.fromFile(self.GSHHS_f_vector_filename)
            gradationRaster_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
            gradationRaster_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shoreline.setRasterResolution(800,800)
            gradationRaster_shoreline.setGradationParameters(150.0,15000.0,1.0)
            gradationRaster_shoreline.calculateLinearGradation()
            gradationRaster_shoreline.writeNetCDF(self.GSHHS_f_gradation_raster_filename)
            #Create raster for mesh gradation towards Shetlands shorelines.
            shetlands_shorelines = qmesh3.vector.Shapes()
            shetlands_shorelines.fromFile(self.shetlands_shoreline_vector_filename)
            gradationRaster_shetlands_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shetlands_shoreline.setShapes(shetlands_shorelines)
            gradationRaster_shetlands_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shetlands_shoreline.setRasterResolution(800,800)
            gradationRaster_shetlands_shoreline.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_shetlands_shoreline.calculateLinearGradation()
            gradationRaster_shetlands_shoreline.writeNetCDF(
                self.shetlands_shoreline_gradation_raster)
            #Create raster for mesh gradation towards 0m gebco contour on the Scottish mainland.
            scotlald_shoreline_GEBCO08 = qmesh3.vector.Shapes()
            scotlald_shoreline_GEBCO08.fromFile(self.scotlald_shoreline_GEBCO08_vector_filename)
            gradationRaster_GEBCO08_0mContour = qmesh3.raster.gradationToShapes()
            gradationRaster_GEBCO08_0mContour.setShapes(scotlald_shoreline_GEBCO08)
            gradationRaster_GEBCO08_0mContour.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GEBCO08_0mContour.setRasterResolution(800,800)
            gradationRaster_GEBCO08_0mContour.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GEBCO08_0mContour.calculateLinearGradation()
            gradationRaster_GEBCO08_0mContour.writeNetCDF(self.gradation_GEBCO08_0mContour_filename)
            #Create raster for mesh gradation towards GSHHS high-resolution lines on the Scottish mainland coast
            GSHHS_h_L1_lines = qmesh3.vector.Shapes()
            GSHHS_h_L1_lines.fromFile(self.GSHHS_h_vector_filename)
            gradationRaster_GSHHS_h_L1_lines = qmesh3.raster.gradationToShapes()
            gradationRaster_GSHHS_h_L1_lines.setShapes(GSHHS_h_L1_lines)
            gradationRaster_GSHHS_h_L1_lines.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GSHHS_h_L1_lines.setRasterResolution(800,800)
            gradationRaster_GSHHS_h_L1_lines.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GSHHS_h_L1_lines.calculateLinearGradation()
            gradationRaster_GSHHS_h_L1_lines.writeNetCDF(self.GSHHS_h_gradation_raster_filename)
            #Create raster for mesh gradation towards Inner Sound plot
            gradationRaster_innerSound_plot = qmesh3.raster.gradationToShapes()
            gradationRaster_innerSound_plot.setShapes(innerSound_plot_polygon)
            gradationRaster_innerSound_plot.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_innerSound_plot.setRasterResolution(800,800)
            gradationRaster_innerSound_plot.setGradationParameters(10.0,15000.0,1.0,0.01)
            gradationRaster_innerSound_plot.calculateLinearGradation()
            gradationRaster_innerSound_plot.writeNetCDF(self.gradation_to_InnerSoundPlot_filename)
            #Calculate overall mesh-metric raster
            meshMetricRaster = qmesh3.raster.minimumRaster([gradationRaster_shoreline, \
                                                           gradationRaster_shetlands_shoreline, \
                                                           gradationRaster_GEBCO08_0mContour, \
                                                           gradationRaster_GSHHS_h_L1_lines, \
                                                           gradationRaster_innerSound_plot])
            meshMetricRaster.writeNetCDF(self.meshMetric_raster_filename)
            #Generate mesh
            domain = qmesh3.mesh.Domain()
            domain.setGeometry(loopShapes, polygonShapes)
            domain.setMeshMetricField(meshMetricRaster)
            domain.setTargetCoordRefSystem('PCC', fldFillValue=1000.0)
            domain.gmsh(geoFilename = self.PCC_geo_filename, \
                        fldFilename = self.PCC_fld_filename, \
                        mshFilename = self.PCC_msh_filename)
        except AssertionError:
            self.assertTrue(False)

    def test_OrkneyShetlandIsles_PCC_tidalSites(self):
        '''Todo: add docstring '''
        try:
            #Try importing qmesh. If an ImportError exception is thrown, then
            # the test might be run in a development setting, and see if qmesh
            # can be found at the project root.
            try:
                import qmesh3
            except ImportError:
                sys.path.append(self.project_root_path)
                import qmesh3
            #Initialise qmesh
            qmesh3.initialise()
            #Read-in the shapefile describing the domain boundaries, and creating a gmsh file.
            boundaries = qmesh3.vector.Shapes()
            boundaries.fromFile(self.all_boundaries_vector_filename)
            loopShapes = qmesh3.vector.identifyLoops(boundaries, isGlobal=False,
                                                    defaultPhysID=1000, fixOpenLoops=True)
            polygonShapes = qmesh3.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=50000)
            #Create loops and polygons for Inner Sound
            innerSound_plot_lines = qmesh3.vector.Shapes()
            innerSound_plot_lines.fromFile(self.innerSound_plot_vector_filename)
            innerSound_plot_loops = qmesh3.vector.identifyLoops(innerSound_plot_lines,
                                                               fixOpenLoops=True)
            innerSound_plot_polygon = qmesh3.vector.identifyPolygons(innerSound_plot_loops, 
                                                                    meshedAreaPhysID = 2)
            #Insert tidal plots and turbines into domain.
            domainLines, domainPolygons = \
               qmesh3.vector.insertRegions(
                  loopShapes, polygonShapes,\
                  innerSound_plot_loops, innerSound_plot_polygon)
            #Create raster for mesh gradation towards full-resolution shorelines.
            GSHHS_fine_boundaries = qmesh3.vector.Shapes()
            GSHHS_fine_boundaries.fromFile(self.GSHHS_f_vector_filename)
            gradationRaster_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
            gradationRaster_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shoreline.setRasterResolution(800,800)
            gradationRaster_shoreline.setGradationParameters(150.0,15000.0,1.0)
            gradationRaster_shoreline.calculateLinearGradation()
            gradationRaster_shoreline.writeNetCDF(self.GSHHS_f_gradation_raster_filename)
            #Create raster for mesh gradation towards Shetlands shorelines.
            shetlands_shorelines = qmesh3.vector.Shapes()
            shetlands_shorelines.fromFile(self.shetlands_shoreline_vector_filename)
            gradationRaster_shetlands_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shetlands_shoreline.setShapes(shetlands_shorelines)
            gradationRaster_shetlands_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shetlands_shoreline.setRasterResolution(800,800)
            gradationRaster_shetlands_shoreline.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_shetlands_shoreline.calculateLinearGradation()
            gradationRaster_shetlands_shoreline.writeNetCDF(
                self.shetlands_shoreline_gradation_raster)
            #Create raster for mesh gradation towards 0m gebco contour on the Scottish mainland.
            scotlald_shoreline_GEBCO08 = qmesh3.vector.Shapes()
            scotlald_shoreline_GEBCO08.fromFile(self.scotlald_shoreline_GEBCO08_vector_filename)
            gradationRaster_GEBCO08_0mContour = qmesh3.raster.gradationToShapes()
            gradationRaster_GEBCO08_0mContour.setShapes(scotlald_shoreline_GEBCO08)
            gradationRaster_GEBCO08_0mContour.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GEBCO08_0mContour.setRasterResolution(800,800)
            gradationRaster_GEBCO08_0mContour.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GEBCO08_0mContour.calculateLinearGradation()
            gradationRaster_GEBCO08_0mContour.writeNetCDF(self.gradation_GEBCO08_0mContour_filename)
            #Create raster for mesh gradation towards GSHHS high-resolution lines on the Scottish mainland coast
            GSHHS_h_L1_lines = qmesh3.vector.Shapes()
            GSHHS_h_L1_lines.fromFile(self.GSHHS_h_vector_filename)
            gradationRaster_GSHHS_h_L1_lines = qmesh3.raster.gradationToShapes()
            gradationRaster_GSHHS_h_L1_lines.setShapes(GSHHS_h_L1_lines)
            gradationRaster_GSHHS_h_L1_lines.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GSHHS_h_L1_lines.setRasterResolution(800,800)
            gradationRaster_GSHHS_h_L1_lines.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GSHHS_h_L1_lines.calculateLinearGradation()
            gradationRaster_GSHHS_h_L1_lines.writeNetCDF(self.GSHHS_h_gradation_raster_filename)
            #Create raster for mesh gradation towards Inner Sound plot
            gradationRaster_innerSound_plot = qmesh3.raster.gradationToShapes()
            gradationRaster_innerSound_plot.setShapes(innerSound_plot_polygon)
            gradationRaster_innerSound_plot.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_innerSound_plot.setRasterResolution(800,800)
            gradationRaster_innerSound_plot.setGradationParameters(10.0,15000.0,1.0,0.01)
            gradationRaster_innerSound_plot.calculateLinearGradation()
            gradationRaster_innerSound_plot.writeNetCDF(self.gradation_to_InnerSoundPlot_filename)
            #Calculate overall mesh-metric raster
            meshMetricRaster = qmesh3.raster.minimumRaster([gradationRaster_shoreline, \
                                                           gradationRaster_shetlands_shoreline, \
                                                           gradationRaster_GEBCO08_0mContour, \
                                                           gradationRaster_GSHHS_h_L1_lines, \
                                                           gradationRaster_innerSound_plot])
            meshMetricRaster.writeNetCDF(self.meshMetric_raster_filename)
            #Generate mesh
            domain = qmesh3.mesh.Domain()
            domain.setGeometry(domainLines, domainPolygons)
            domain.setMeshMetricField(meshMetricRaster)
            domain.setTargetCoordRefSystem('PCC', fldFillValue=1000.0)
            domain.gmsh(geoFilename = self.PCC_inner_sound_plot_geo_filename, \
                        fldFilename = self.PCC_inner_sound_plot_fld_filename, \
                        mshFilename = self.PCC_inner_sound_plot_msh_filename)
        except AssertionError:
            self.assertTrue(False)

    def test_OrkneyShetlandIsles_PCC_InnerSoundTurbines(self):
        '''Todo: add docstring '''
        try:
            #Try importing qmesh. If an ImportError exception is thrown, then
            # the test might be run in a development setting, and see if qmesh
            # can be found at the project root.
            try:
                import qmesh3
            except ImportError:
                sys.path.append(self.project_root_path)
                import qmesh3
            import qgis.core
            import numpy
            qmesh3.initialise()
            #Read-in the shapefile describing the domain boundaries, and creating a gmsh file.
            boundaries = qmesh3.vector.Shapes()
            boundaries.fromFile(self.all_boundaries_vector_filename)
            loopShapes = qmesh3.vector.identifyLoops(boundaries, isGlobal=False,
                                                    defaultPhysID=1000, fixOpenLoops=True)
            polygonShapes = qmesh3.vector.identifyPolygons(loopShapes, smallestNotMeshedArea=50000)
            #Create loops and polygons for Inner Sound
            innerSound_plot_lines = qmesh3.vector.Shapes()
            innerSound_plot_lines.fromFile(self.innerSound_plot_vector_filename)
            innerSound_plot_loops = qmesh3.vector.identifyLoops(innerSound_plot_lines,
                                                               fixOpenLoops=True)
            innerSound_plot_polygon = qmesh3.vector.identifyPolygons(innerSound_plot_loops, 
                                                                    meshedAreaPhysID = 2)
            #Using a random-number generator to distribute points inside
            # the Inner Sound plot, with the additional constrains that
            # they must be 200m apart. Locate 180 such points. Then create
            # a cicle for each point, representative of a tidal turbine.
            #Seed the random-number generator so as to obtain same set of
            # points every time the test is run. 
            numpy.random.seed([1,0])
            #Calculate mapping
            # parameters so as to map the interval [0,1] to an area
            # enclosing the tidal plot. The numbers below are corner
            # coordinates from a rectangle -in UTM30- that ecloses
            # the tidal plot
            ksi_param_1 = 494800 - 490400
            ksi_param_2 = 490400
            eta_param_1 = 6503800 - 6501600
            eta_param_2 = 6501600
            turbinePoints = []
            #Change the coordinate reference system of the Inner Sound plot
            # from lon-lat into UTM30 .
            innerSound_plot_polygon.changeCoordRefSystem('EPSG:32630')
            polygon = innerSound_plot_polygon.getFeatures()[0]
            while len(turbinePoints) < 180:
                #Get new "arbitrary" point, carry out mapping using
                # the parameters calculated above.
                newPointX = numpy.random.random()*ksi_param_1 + ksi_param_2
                newPointY = numpy.random.random()*eta_param_1 + eta_param_2
                #Construct a large circle around that point, 15m in radius
                # Then check if it is included into the tidal plot. The larger
                # radius should eliminate the possibility of the centre-point
                # of the circle lying inside the polygon, but too close to the
                # plot boundary such that a part of the tidal-turbine is
                # outside the plot.
                turbineBigCircle = qmesh3.vector.Circle(newPointX, newPointY, 15.0, 20, 'EPSG:32630')
                turbineBigPolygon = turbineBigCircle.asQgsPolygonFeature()
                if polygon.geometry().contains(turbineBigPolygon.geometry()):
                    if len(turbinePoints) == 0:
                        turbinePoints.append((newPointX, newPointY))
                    else:
                        isNewPointAdmissible = True
                        for otherPoint in turbinePoints:
                            distance = numpy.sqrt( (otherPoint[0] - newPointX)**2 + (otherPoint[1] - newPointY)**2)
                            if distance < 100.0:
                                isNewPointAdmissible = False
                                break
                        if isNewPointAdmissible:
                            turbinePoints.append((newPointX, newPointY))
            turbines = qmesh3.vector.Shapes()
            turbines.setCoordRefSystemFromString('EPSG:32630')
            turbines.setShapeType(qgis.core.QgsWkbTypes.LineString)
            for turbinePoint in turbinePoints:
                turbineCircle = qmesh3.vector.Circle(turbinePoint[0], turbinePoint[1], 10.0, 20, 'EPSG:32630')
                turbines.addFeature(turbineCircle.asQgsFeature())
            turbines.writeFile(self.turbine_lines_filename)
            turbines.changeCoordRefSystem('EPSG:4326')
            turbineLoops = qmesh3.vector.identifyLoops(turbines,
                      fixOpenLoops=True)
            turbinePolygons = qmesh3.vector.identifyPolygons(turbineLoops,
                      meshedAreaPhysID = 3)
            #Reconstruct innerSound_plot_polygon object, so as to obtain
            # it back to EPSG:4326, recall we changed its CRS to EPSG:32630
            # While we could invoke the changeCoordRefSystem() method on it
            # we will not obtain the same point coordinates as before, leading
            # to problems later on.
            innerSound_plot_polygon = qmesh3.vector.identifyPolygons(innerSound_plot_loops, 
                                                                    meshedAreaPhysID = 2)
            #Insert tidal plots and turbines into domain.
            domainLines, domainPolygons = \
               qmesh3.vector.insertRegions(
                  loopShapes, polygonShapes,\
                  innerSound_plot_loops, innerSound_plot_polygon)
            domainLines, domainPolygons = \
               qmesh3.vector.insertRegions(
                  domainLines, domainPolygons,\
                  turbineLoops, turbinePolygons)
            domainLines.writeFile(self.domain_lines_filename)
            domainPolygons.writeFile(self.domain_polygons_filename)
            #Create raster for mesh gradation towards full-resolution shorelines.
            GSHHS_fine_boundaries = qmesh3.vector.Shapes()
            GSHHS_fine_boundaries.fromFile(self.GSHHS_f_vector_filename)
            gradationRaster_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shoreline.setShapes(GSHHS_fine_boundaries)
            gradationRaster_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shoreline.setRasterResolution(800,800)
            gradationRaster_shoreline.setGradationParameters(150.0,15000.0,1.0)
            gradationRaster_shoreline.calculateLinearGradation()
            gradationRaster_shoreline.writeNetCDF(self.GSHHS_f_gradation_raster_filename)
            #Create raster for mesh gradation towards Shetlands shorelines.
            shetlands_shorelines = qmesh3.vector.Shapes()
            shetlands_shorelines.fromFile(self.shetlands_shoreline_vector_filename)
            gradationRaster_shetlands_shoreline = qmesh3.raster.gradationToShapes()
            gradationRaster_shetlands_shoreline.setShapes(shetlands_shorelines)
            gradationRaster_shetlands_shoreline.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_shetlands_shoreline.setRasterResolution(800,800)
            gradationRaster_shetlands_shoreline.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_shetlands_shoreline.calculateLinearGradation()
            gradationRaster_shetlands_shoreline.writeNetCDF(
                self.shetlands_shoreline_gradation_raster)
            #Create raster for mesh gradation towards 0m gebco contour on the Scottish mainland.
            scotlald_shoreline_GEBCO08 = qmesh3.vector.Shapes()
            scotlald_shoreline_GEBCO08.fromFile(self.scotlald_shoreline_GEBCO08_vector_filename)
            gradationRaster_GEBCO08_0mContour = qmesh3.raster.gradationToShapes()
            gradationRaster_GEBCO08_0mContour.setShapes(scotlald_shoreline_GEBCO08)
            gradationRaster_GEBCO08_0mContour.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GEBCO08_0mContour.setRasterResolution(800,800)
            gradationRaster_GEBCO08_0mContour.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GEBCO08_0mContour.calculateLinearGradation()
            gradationRaster_GEBCO08_0mContour.writeNetCDF(self.gradation_GEBCO08_0mContour_filename)
            #Create raster for mesh gradation towards GSHHS high-resolution lines on the Scottish mainland coast
            GSHHS_h_L1_lines = qmesh3.vector.Shapes()
            GSHHS_h_L1_lines.fromFile(self.GSHHS_h_vector_filename)
            gradationRaster_GSHHS_h_L1_lines = qmesh3.raster.gradationToShapes()
            gradationRaster_GSHHS_h_L1_lines.setShapes(GSHHS_h_L1_lines)
            gradationRaster_GSHHS_h_L1_lines.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_GSHHS_h_L1_lines.setRasterResolution(800,800)
            gradationRaster_GSHHS_h_L1_lines.setGradationParameters(1500.0,15000.0,0.5)
            gradationRaster_GSHHS_h_L1_lines.calculateLinearGradation()
            gradationRaster_GSHHS_h_L1_lines.writeNetCDF(self.GSHHS_h_gradation_raster_filename)
            #Create raster for mesh gradation towards Inner Sound plot
            gradationRaster_innerSound_plot = qmesh3.raster.gradationToShapes()
            gradationRaster_innerSound_plot.setShapes(innerSound_plot_polygon)
            gradationRaster_innerSound_plot.setRasterBounds(-6.0, 1.0, 57.0, 62.0)
            gradationRaster_innerSound_plot.setRasterResolution(800,800)
            gradationRaster_innerSound_plot.setGradationParameters(10.0,15000.0,1.0,0.01)
            gradationRaster_innerSound_plot.calculateLinearGradation()
            gradationRaster_innerSound_plot.writeNetCDF(self.gradation_to_InnerSoundPlot_filename)
            #Calculate overall mesh-metric raster
            meshMetricRaster = qmesh3.raster.minimumRaster([gradationRaster_shoreline, \
                                                           gradationRaster_shetlands_shoreline, \
                                                           gradationRaster_GEBCO08_0mContour, \
                                                           gradationRaster_GSHHS_h_L1_lines, \
                                                           gradationRaster_innerSound_plot])
            meshMetricRaster.writeNetCDF(self.meshMetric_raster_filename)
            #Generate mesh
            domain = qmesh3.mesh.Domain()
            domain.setGeometry(domainLines, domainPolygons)
            domain.setMeshMetricField(meshMetricRaster)
            domain.setTargetCoordRefSystem('PCC', fldFillValue=1000.0)
            domain.gmsh(geoFilename = self.PCC_inner_sound_turbines_geo_filename, \
                        fldFilename = self.PCC_inner_sound_turbines_fld_filename, \
                        mshFilename = self.PCC_inner_sound_turbines_msh_filename)
            #Convert mesh into shapefile-set. The mesh vertices are re-projected, which produces
            # a warning. It is safe to ignore the warning, as we do the reprojection for
            # visualisation purposes alone.
            mesh = qmesh3.mesh.Mesh()
            mesh.readGmsh(self.PCC_inner_sound_turbines_msh_filename, 'PCC')
            mesh.reProjectVertices('EPSG:4326')
            mesh.writeShapefile(self.PCC_inner_sound_turbines_shp_filename)
        except AssertionError:
            self.assert_(False)
 
if __name__ == '__main__':
    unittest.main()
