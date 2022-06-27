import numpy as np
import qgis.core
import PyQt4.QtCore
import loxodrome

outputFilename = PyQt4.QtCore.QString('blendedLoxodromes.shp')

qgis.core.QgsApplication.setPrefixPath('/usr', True)
qgis.core.QgsApplication.initQgis()

fields = {'ID':qgis.core.QgsField("ID", PyQt4.QtCore.QVariant.Int),
          'Descr': qgis.core.QgsField("Descrip", PyQt4.QtCore.QVariant.String) }

writer = qgis.core.QgsVectorFileWriter(outputFilename, "CP1250", fields, qgis.core.QGis.WKBLineString, None, "ESRI Shapefile")
#writer = qgis.core.QgsVectorFileWriter(outputFilename, "CP1250", fields, qgis.core.QGis.WKBPoint, None, "ESRI Shapefile")
if writer.hasError() != qgis.core.QgsVectorFileWriter.NoError:
  print("Error when creating shapefile: ", writer.hasError())

#Create the loxodromes.
ldrome_1 = loxodrome.loxodromicLine([51.0,10.5], 130.0)
ldrome_1_segment = ldrome_1.segment(100,0, ldrome_1.point, [54.4, np.nan])
ldrome_2 = loxodrome.loxodromicLine([54.4,17.5], 150.0)
ldrome_2_segment = ldrome_2.segment(100,0, ldrome_2.point, [np.nan, 10.5])
ldrome_2_segment.reverse()
#Blend the loxodromes, to create the open boundary.
segment=[]
ldrome_1_weight = np.linspace(1,0,len(ldrome_1_segment))
ldrome_2_weight = np.linspace(0,1,len(ldrome_2_segment))
for index in range(len(ldrome_1_segment)):
    x = ldrome_1_segment[index][0]*ldrome_1_weight[index] + ldrome_2_segment[index][0]*ldrome_2_weight[index]
    y = ldrome_1_segment[index][1]*ldrome_1_weight[index] + ldrome_2_segment[index][1]*ldrome_2_weight[index]
    segment.append([x,y])

#Convert blended loxodromes into a list if QgsPoint, ready for witing to shapefile.
segment_Qgs=[]
for point in segment:
    segment_Qgs.append(qgis.core.QgsPoint(point[0],point[1]))

#Create an empty feature, set the geometry using the loxodrome points and write. 
feat = qgis.core.QgsFeature()
feat.setGeometry(qgis.core.QgsGeometry.fromPolyline(segment_Qgs))
feat.setAttributeMap({'ID':1000,'Descr':'Blended loxodromes'})
writer.addFeature(feat)

# delete the writer to flush features to disk (optional)
del writer
