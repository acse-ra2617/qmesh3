import numpy as np
import qmesh

qmesh.initialise()

#Create the loxodromes at the Western side of the domain
ldrome_1 = qmesh.vector.loxodromicLine((-5.0,58.6), -10.0, 100, 0, (np.nan, 60.2), "EPSG:4326")
ldrome_2 = qmesh.vector.loxodromicLine((-3.43,60.46), -135.0, 100, 0, (-5.0, np.nan), "EPSG:4326")
#Blend the loxodromes, to create the open boundary.
openBoundary = qmesh.vector.blendedLoxodromes(ldrome_1,ldrome_2,100)
#Write to file
openBoundary.asShapes().writeFile('blendedLoxodromesWest.shp')

#Create the loxodromes at the Eastern side of the domain
ldrome_1 = qmesh.vector.loxodromicLine((-1.83,57.6), 45.0, 100, 0, (np.nan, 61.94), "EPSG:4326")
ldrome_2 = qmesh.vector.loxodromicLine((0.0,61.94), 90.0, 100, 0, (2.0, np.nan), "EPSG:4326")
#Blend the loxodromes, to create the open boundary.
openBoundary = qmesh.vector.blendedLoxodromes(ldrome_1,ldrome_2,100)
#Write to file
openBoundary.asShapes().writeFile('blendedLoxodromesEast.shp')
