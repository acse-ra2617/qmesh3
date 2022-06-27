===========================
Tutorials
===========================

The Orkney and Shetland Islands
-------------------------------

A bit about the region

Introduce the two basic data-structures of GIS and talk about their correspondence to mesh generator structures.

Show one we did earlier: Show and discuss figures of the vector and rasters used for the Orkneys example. Show and discuss figures of the mesh

Introduce QGIS, introduce shapefile and raster representation as layers, introduce canvas and Coordinate Reference Systems. Introduce Gmsh. Give a few links and citations for both fostware.

The rest of this tutorial delineates creation of the mesh from the vector data. We are assuming a linux. Text in ``literal`` prepended by ``$`` below signify commands to be typed into a linux terminal.


Building the geometry: shorelines and open boundaries.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For this part of the tutorial we will fetch an existing qgis project containing a near complete domain geometry description. The reader will be instructed how to prepare the geometry for meshing.

Download a data-set containing the qgis project file, bathymetry NetCDF and shapefiles from `figshare <http://figshare.com/s/f4d8427eb9cc11e4b74f06ec4bbcf141>`_: Clicking on the "Download all" button will download a zip file of the data-set. Move the zip file to a loaction of your choosing, preferably a directory, and unizip 

``$ unzip 1314221.zip`` 

Open the Qgis project file with Qgis 

``$ qgis OrkneyShetlandIsles.qgs &``

You should see something similar to figure:qgisWindowOrkneyShetlands_

  .. _figure:qgisWindowOrkneyShetlands
  .. figure :: figures/qgisWindowOrkneyShetlands.png
     :align: center
     :scale: 75 %
     :figclass: align-center
  The Qgis Window, showing the various shapefiles and bathymetry raster of the Orkney and Shetlands example.

The coloured lines in figure:qgisWindowOrkneyShetlands_ will be used to assemble the geometry definition of the domain. Notice that data originates from a different sources, as shown in the left panel of the qgis window. The black lines originate from full-resolution GSHHG data while the purple lines originate from the "high" resolution GSHHG data (add GSHHG citation). The red and blue lines were extracted from the GEBCO08 2012 bathymetry. The red line is a zero-meter elevation contour and the blue line is a contour at a depth of 300m. The yellow lines are arbitrarily drawn, to close the domain. The bathymetry raster is an excerpt from the GEBCO08 2012 bathymetry (add gebco citation)

As discussed above the lines are classified as "vector" and are encoded as shapefile-sets, while the bathymetry is a "raster" and here is encoded as a NetCDF file. Qgis represents each shapefile set as a separate "vector layer" and raster data is represented as separate "raster layers". Note that Qgis has functionality that allows contours to be extracted from raster data. It is left as an exercise for the reader to discover how to do such extractions (hint: click on Raster -> Extraction -> Contour). This contour extraction functionality was used to construct the red and blue lines. In this way the Qgis user can combine information from a huge variety of sources to create geometry definitions for complex realistic ocean domains.

The domain is loosely defined. If you notice the ends of the yellow lines they do not exactly meet the end-points of other lines. Yet, Qmesh needs all lines to meet so that they precisely define the domain. Thankfully, such vector editing is simple with Qgis. Lets go through the steps

First collect a copy all lines onto a single layer:

   #. Create a new vector layer: Click on the "Layer" menu on the toolbar and select "Create Layer", then "New Shapefile Layer..." from the drop-down menus. A similar dialogue window to that in figure:newShapefileLayerDialogue_ will appear

   #. On the top row of the dialogue window select "line". Click on next, you will be prompted to save the shapefile. Name your shapefile "domainBoundaries" and make sure you save it in the same directory as the other files of this tutorial, as shown in figure:storeDomainBoundariesShapefile_. When you click save the shapefile is created, and a vector layer named "domainBoundaries" appears in the Qgis window.

   #. Select the "shelf_break_300mContour" layer on the layers side panel. Click on the "Select Features" tool on the toolbar, highlighted in figure:qgisWindowOrkneyShetlands_selectorTool_ . Notice the downwards pointing triangle next to the selector tool, clicking on that will allow you to use different variants of the tool, as shown in figure figure:qgisWindowOrkneyShetlands_selectorTool_. Choose the "select Features by Freehand" tool.

   #. Click and drag across the canvas, the selector tool will outline an area. Any features on the current layer with any point in the highlighted area will be selected.


  .. _figure:newShapefileLayerDialogue
  .. figure :: figures/newShapefileLayerDialogue.png
     :align: center
     :scale: 75 %
     :figclass: align-center
  The Qgis new shapefile dialogue window.

  .. _figure:storeDomainBoundariesShapefile
  .. figure :: figures/storeDomainBoundariesShapefile.png
     :align: center
     :scale: 75 %
     :figclass: align-center
  The Qgis "Save layer as..." dialogue window.

  .. _figure:qgisWindowOrkneyShetlands_selectorTool
  .. figure :: figures/qgisWindowOrkneyShetlands_selectorTool.png
     :align: center
     :scale: 75 %
     :figclass: align-center
  The Qgis main window, the "shelf_break_300mContour" layer and the "select Features by Freehand" tool have been selected.


Defining the mesh size
^^^^^^^^^^^^^^^^^^^^^^

Discuss ideal mesh spacing for this example.

Walk reader through mesh metric creation.

Meshing
^^^^^^^

Get user to use the qmesh GUI to mesh the project, possibly in more than one CRS.

