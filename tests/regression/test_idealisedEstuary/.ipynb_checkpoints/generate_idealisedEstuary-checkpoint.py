import numpy as np
import qgis.core
import PyQt5.QtCore
import loxodrome
import sys
sys.path.insert(0,"../../../")
from qmesh3 import vector

outputLineFilename = 'idealisedEstuary_lines.shp'
outputPolygonFilename = 'idealisedEstuary_polygon.shp'

qgis.core.QgsApplication.setPrefixPath('/usr', True)
qgis.core.QgsApplication.initQgis()

IDfield = qgis.core.QgsField("ID", PyQt5.QtCore.QVariant.Int)
descriptionField = qgis.core.QgsField("Descript", PyQt5.QtCore.QVariant.String)
physicalIDfield = qgis.core.QgsField("PhysID", PyQt5.QtCore.QVariant.Int)
fields = qgis.core.QgsFields()
fields.append(IDfield)
fields.append(descriptionField)
fields.append(physicalIDfield)

coordinateReferenceSystem = qgis.core.QgsCoordinateReferenceSystem()
coordinateReferenceSystem.createFromString("EPSG:4326")

writer = qgis.core.QgsVectorFileWriter(outputLineFilename, "CP1250", fields, qgis.core.QgsWkbTypes.LineString, coordinateReferenceSystem, "ESRI Shapefile")
if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
  print("Error when creating shapefile: ", writer.hasError())

#Create the boundaries.
#Starting with the West-side boundary, create segment between
# (-10,-5) and (-15,0) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([-10.0, -15.0], 0.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [np.nan, 0.0])
ldrome_2 = loxodrome.loxodromicLine([-15.0, 0.0], 160.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [np.nan, -15.0])
ldrome_2_segment.reverse()
segment=[]
ldrome_1_weight = (1. + np.cos(np.linspace(0,np.pi,len(ldrome_1_segment))))/2.
ldrome_2_weight = (1. + np.cos(np.linspace(np.pi,0,len(ldrome_2_segment))))/2.
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature1 = qgis.core.QgsFeature()
feature1.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature1.setFields(fields)
feature1.setAttribute('ID',1)
feature1.setAttribute('Descript','Lower West boundary segment')
feature1.setAttribute('PhysID',1000)
writer.addFeature(feature1)

#On the West-side boundary, create segment between
# (-15,0) and (-25,30) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([-15.0, 0.0], -20.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [np.nan, 30.0])
ldrome_2 = loxodrome.loxodromicLine([-25.0, 30.0], 130.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [np.nan, 0.0])
ldrome_2_segment.reverse()
segment=[]
ldrome_1_weight = (1. + np.cos(np.linspace(0,np.pi,len(ldrome_1_segment))))/2.
ldrome_2_weight = (1. + np.cos(np.linspace(np.pi,0,len(ldrome_2_segment))))/2.
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature2 = qgis.core.QgsFeature()
feature2.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature2.setFields(fields)
feature2.setAttribute('ID',2)
feature2.setAttribute('Descript','Upper West boundary segment')
feature2.setAttribute('PhysID',1000)
writer.addFeature(feature2)

#On the East-side boundary, create segment between
# (10,-15) and (15,0) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([10.0, -15.0], 0.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [np.nan, 0.0])
ldrome_2 = loxodrome.loxodromicLine([15.0, 0.0], -160.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [np.nan, -15.0])
ldrome_2_segment.reverse()
segment=[]
ldrome_1_weight = (1. + np.cos(np.linspace(0,np.pi,len(ldrome_1_segment))))/2.
ldrome_2_weight = (1. + np.cos(np.linspace(np.pi,0,len(ldrome_2_segment))))/2.
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature3 = qgis.core.QgsFeature()
feature3.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature3.setFields(fields)
feature3.setAttribute('ID',3)
feature3.setAttribute('Descript','Lower East boundary segment')
feature3.setAttribute('PhysID',1000)
writer.addFeature(feature3)

#On the East-side boundary, create segment between
# (15,0) and (25,30) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([15.0, 0.0], 20.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [np.nan, 30.0])
ldrome_2 = loxodrome.loxodromicLine([25.0, 30.0], -130.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [np.nan, 0.0])
ldrome_2_segment.reverse()
segment=[]
ldrome_1_weight = (1. + np.cos(np.linspace(0,np.pi,len(ldrome_1_segment))))/2.
ldrome_2_weight = (1. + np.cos(np.linspace(np.pi,0,len(ldrome_2_segment))))/2.
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature4 = qgis.core.QgsFeature()
feature4.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature4.setFields(fields)
feature4.setAttribute('ID',4)
feature4.setAttribute('Descript','Upper East boundary segment')
feature4.setAttribute('PhysID',1000)
writer.addFeature(feature4)

#On the North-side boundary, create segment between
# (-25,30) and (25,30) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([-25.0, 30.0], 60.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [25.0, np.nan])
ldrome_2 = loxodrome.loxodromicLine([25.0, 30.0], -60.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [-25.0, np.nan])
ldrome_2_segment.reverse()
segment=[]
ldrome_1_weight = (1. + np.cos(np.linspace(0,np.pi,len(ldrome_1_segment))))/2.
ldrome_2_weight = (1. + np.cos(np.linspace(np.pi,0,len(ldrome_2_segment))))/2.
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature5 = qgis.core.QgsFeature()
feature5.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature5.setFields(fields)
feature5.setAttribute('ID',5)
feature5.setAttribute('Descript','North boundary segment')
feature5.setAttribute('PhysID',2000)
writer.addFeature(feature5)

#On the South-side boundary, create segment between
# (-10,-15) and (10,-15) as two blended loxodromes.
ldrome_1 = loxodrome.loxodromicLine([-10.0, -15.0], 90.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [10.0, -15])
#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in ldrome_1_segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
#Create an empty feature, set the geometry using the loxodrome points and write. 
feature6 = qgis.core.QgsFeature()
feature6.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature6.setFields(fields)
feature6.setAttribute('ID',6)
feature6.setAttribute('Descript','South boundary segment')
feature6.setAttribute('PhysID',3000)
writer.addFeature(feature6)

#Create a small circular island, composed of two lines
islandCenter = [-5,15]
islandRadius = 2.5
npoints = 50
westSegment=[]
deltaAngle = np.radians(180)/(npoints-1)
for index in range(npoints-1):
    angle = np.radians(-90) - deltaAngle*index
    x = np.cos(angle)*islandRadius + islandCenter[0]
    y = np.sin(angle)*islandRadius + islandCenter[1]
    westSegment.append([x,y])
westSegment.append([islandCenter[0], islandCenter[1]+islandRadius])
segment_Qgs=[]
for point in westSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature7 = qgis.core.QgsFeature()
feature7.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature7.setFields(fields)
feature7.setAttribute('ID',7)
feature7.setAttribute('Descript','Island west boundary segment')
feature7.setAttribute('PhysID',1000)
writer.addFeature(feature7)
eastSegment=[]
for index in range(npoints-1):
    angle = np.radians(-90) + deltaAngle*index
    x = np.cos(angle)*islandRadius + islandCenter[0]
    y = np.sin(angle)*islandRadius + islandCenter[1]
    eastSegment.append([x,y])
eastSegment.append([islandCenter[0], islandCenter[1]+islandRadius])
segment_Qgs=[]
for point in eastSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature8 = qgis.core.QgsFeature()
feature8.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature8.setFields(fields)
feature8.setAttribute('ID',8)
feature8.setAttribute('Descript','Island east boundary segment')
feature8.setAttribute('PhysID',1000)
writer.addFeature(feature8)

#Create another small circular island, composed of two lines
islandCenter = [5,10]
islandRadius = 5.0
npoints = 50
westSegment=[]
deltaAngle = np.radians(180)/(npoints-1)
for index in range(npoints-1):
    angle = np.radians(90) + deltaAngle*index
    x = np.cos(angle)*islandRadius + islandCenter[0]
    y = np.sin(angle)*islandRadius + islandCenter[1]
    westSegment.append([x,y])
westSegment.append([islandCenter[0], islandCenter[1]-islandRadius])
segment_Qgs=[]
for point in westSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature9 = qgis.core.QgsFeature()
feature9.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature9.setFields(fields)
feature9.setAttribute('ID',9)
feature9.setAttribute('Descript','Island west boundary segment')
feature9.setAttribute('PhysID',1000)
writer.addFeature(feature9)
eastSegment=[]
for index in range(npoints-1):
    angle = np.radians(90) - deltaAngle*index
    x = np.cos(angle)*islandRadius + islandCenter[0]
    y = np.sin(angle)*islandRadius + islandCenter[1]
    eastSegment.append([x,y])
eastSegment.append([islandCenter[0], islandCenter[1]-islandRadius])
segment_Qgs=[]
for point in eastSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature10 = qgis.core.QgsFeature()
feature10.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature10.setFields(fields)
feature10.setAttribute('ID',10)
feature10.setAttribute('Descript','Island east boundary segment')
feature10.setAttribute('PhysID',1000)
writer.addFeature(feature10)

#Create a circular patch, to test the 'region' IDs capability.
patchCenter = [5,30]
patchRadius = 10.0
npoints = 50
westSegment=[]
deltaAngle = np.radians(180)/(npoints-1)
for index in range(npoints-1):
    angle = np.radians(90) + deltaAngle*index
    x = np.cos(angle)*patchRadius + patchCenter[0]
    y = np.sin(angle)*patchRadius + patchCenter[1]
    westSegment.append([x,y])
westSegment.append([patchCenter[0], patchCenter[1]-patchRadius])
segment_Qgs=[]
for point in westSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature11 = qgis.core.QgsFeature()
feature11.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature11.setFields(fields)
feature11.setAttribute('ID',11)
feature11.setAttribute('Descript','Patch west boundary segment')
#No PhysID attribute is set, as we do not want physical id to this line
writer.addFeature(feature11)
eastSegment=[]
for index in range(npoints-1):
    angle = np.radians(90) - deltaAngle*index
    x = np.cos(angle)*patchRadius + patchCenter[0]
    y = np.sin(angle)*patchRadius + patchCenter[1]
    eastSegment.append([x,y])
eastSegment.append([patchCenter[0], patchCenter[1]-patchRadius])
segment_Qgs=[]
for point in eastSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature12 = qgis.core.QgsFeature()
feature12.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature12.setFields(fields)
feature12.setAttribute('ID',12)
feature12.setAttribute('Descript','Patch east boundary segment')
#No PhysID attribute is set, as we do not want physical id to this line
writer.addFeature(feature12)

#Create an arctic lake.
lakeCenter = [0,75]
lakeRadius = 10.0
npoints = 50
westSegment=[]
deltaAngle = np.radians(180)/(npoints-1)
for index in range(npoints-1):
    angle = np.radians(90) + deltaAngle*index
    x = np.cos(angle)*lakeRadius + lakeCenter[0]
    y = np.sin(angle)*lakeRadius + lakeCenter[1]
    westSegment.append([x,y])
westSegment.append([lakeCenter[0], lakeCenter[1]-lakeRadius])
segment_Qgs=[]
for point in westSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature13 = qgis.core.QgsFeature()
feature13.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature13.setFields(fields)
feature13.setAttribute('ID',13)
feature13.setAttribute('Descript','Arctic lake west boundary segment')
feature13.setAttribute('PhysID',1000)
writer.addFeature(feature13)
eastSegment=[]
for index in range(npoints-1):
    angle = np.radians(90) - deltaAngle*index
    x = np.cos(angle)*lakeRadius + lakeCenter[0]
    y = np.sin(angle)*lakeRadius + lakeCenter[1]
    eastSegment.append([x,y])
eastSegment.append([lakeCenter[0], lakeCenter[1]-lakeRadius])
segment_Qgs=[]
for point in eastSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature14 = qgis.core.QgsFeature()
feature14.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature14.setFields(fields)
feature14.setAttribute('ID',14)
feature14.setAttribute('Descript','Arctic lake east boundary segment')
feature14.setAttribute('PhysID',1000)
writer.addFeature(feature14)

#Create an arctic lake.
lakeCenter = [0,-75]
lakeRadius = 10.0
npoints = 50
westSegment=[]
deltaAngle = np.radians(180)/(npoints-1)
for index in range(npoints-1):
    angle = np.radians(90) + deltaAngle*index
    x = np.cos(angle)*lakeRadius + lakeCenter[0]
    y = np.sin(angle)*lakeRadius + lakeCenter[1]
    westSegment.append([x,y])
westSegment.append([lakeCenter[0], lakeCenter[1]-lakeRadius])
segment_Qgs=[]
for point in westSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature15 = qgis.core.QgsFeature()
feature15.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature15.setFields(fields)
feature15.setAttribute('ID',15)
feature15.setAttribute('Descript','Antarctic lake west boundary segment')
feature15.setAttribute('PhysID',1000)
writer.addFeature(feature15)
eastSegment=[]
for index in range(npoints-1):
    angle = np.radians(90) - deltaAngle*index
    x = np.cos(angle)*lakeRadius + lakeCenter[0]
    y = np.sin(angle)*lakeRadius + lakeCenter[1]
    eastSegment.append([x,y])
eastSegment.append([lakeCenter[0], lakeCenter[1]-lakeRadius])
segment_Qgs=[]
for point in eastSegment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))
feature16 = qgis.core.QgsFeature()
feature16.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feature16.setFields(fields)
feature16.setAttribute('ID',16)
feature16.setAttribute('Descript','Antarctic lake east boundary segment')
feature16.setAttribute('PhysID',1000)
writer.addFeature(feature16)

# delete the writer to flush features to disk.
del writer

#Open polygon file
polygonWriter = qgis.core.QgsVectorFileWriter(outputPolygonFilename, "CP1250", fields, qgis.core.QgsWkbTypes.MultiPolygon, coordinateReferenceSystem, "ESRI Shapefile")
if polygonWriter.hasError() != qgis.core.QgsVectorFileWriter.NoError:
  print("Error when creating polygon shapefile: ", writer.hasError())

#Create polygon features

#Create polgon feature for 'estuary'
polygonFeature1 = vector.shapefileTools.lines2polygon([feature1, feature2, feature3, feature4, feature5, feature6], [[feature7, feature8], [feature9, feature10], [feature11, feature12]])
polygonFeature1.setFields(fields)
polygonFeature1.setAttribute('ID',13)
polygonFeature1.setAttribute('Descript','Idealised Estuary Polygon')
polygonFeature1.setAttribute('PhysID',10000)
polygonWriter.addFeature(polygonFeature1)

#Create polygon feature for patch in 'estuary'
polygonFeature2 = vector.shapefileTools.lines2polygon([feature11, feature12], [])
polygonFeature2.setFields(fields)
polygonFeature2.setAttribute('ID',14)
polygonFeature2.setAttribute('Descript','Circular patch')
polygonFeature2.setAttribute('PhysID',20000)
polygonWriter.addFeature(polygonFeature2)

#Create polygon feature for 'Arctic lake'
polygonFeature3 = vector.shapefileTools.lines2polygon([feature13, feature14], [])
polygonFeature3.setFields(fields)
polygonFeature3.setAttribute('ID',15)
polygonFeature3.setAttribute('Descript','Arctic lake')
polygonFeature3.setAttribute('PhysID',30000)
polygonWriter.addFeature(polygonFeature3)

#Create polygon feature for 'Arctic lake'
polygonFeature4 = vector.shapefileTools.lines2polygon([feature15, feature16], [])
polygonFeature4.setFields(fields)
polygonFeature4.setAttribute('ID',14)
polygonFeature4.setAttribute('Descript','Antarctic lake')
polygonFeature4.setAttribute('PhysID',40000)
polygonWriter.addFeature(polygonFeature4)

# delete the writer to flush features to disk.
del polygonWriter
