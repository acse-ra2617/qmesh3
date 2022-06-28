import numpy as np

class loxodromicLine(object):
    def __init__(self, point, trueNorthBearing):
        self.point = point
        self.trueNorthBearing = trueNorthBearing
        self.alpha = np.tan(np.radians(90. - trueNorthBearing))
        self.beta = point[1] - self.alpha*point[0]

    def segment(self, numbPoints, loopArounds, startPoint, endPoint):
        #Check that points, if fully specified, actually lie on
        # the loxodrome.
        if not(np.isnan(startPoint[0])) and\
           not(np.isnan(startPoint[1])):
            if abs(startPoint[1] - self.alpha*startPoint[0] - self.beta) > 1e-12:
                raise Exception
        if not(np.isnan(endPoint[0])) and\
           not(np.isnan(endPoint[1])):
            if abs(endPoint[1] - self.alpha*(endPoint[0] + loopArounds*360.0) - self.beta) > 1e-12:
                raise Exception
        #Calculate coordinates of end-point, if not fully specified.
        if np.isnan(endPoint[0]) and not(np.isnan(endPoint[1])):
            endPoint[0] = ((endPoint[1] - self.beta)/self.alpha) -\
                          loopArounds*360.0
        if not(np.isnan(endPoint[0])) and np.isnan(endPoint[1]):
            endPoint[1] = self.alpha*(endPoint[0] + loopArounds*360.0) +\
                          self.beta
        #Create points along equi-spaced intervals on loxodrome.
        segment = []
        delta_lon = ((endPoint[0] + loopArounds*360.0) - startPoint[0])/(numbPoints - 1)
        #Delta longitude will be 0 if the loxodrome is drawn at a true nort bearing of 0,
        # so lets calculate Delta latitude as well
        delta_lat = (endPoint[1] - startPoint[1])/(numbPoints - 1)
        if delta_lon != 0.0:
            for pointIndex in range(numbPoints):
                lon = startPoint[0] + pointIndex*delta_lon
                lat = self.alpha*lon + self.beta
                lon = np.sign(lon)*np.mod(abs(lon),360.0)
                #if Latitude is greater than 180 or smaller than -180, the loxodrome has reached the pole and must stop.
                if lat >= 180.0:
                    segment.append([lon,180.0])
                    return segment
                if lat <=-180.0:
                    segment.append([lon,-180.0])
                    return segment
                else:
                    segment.append([lon,lat])
        else:
            for pointIndex in range(numbPoints):
                lon = startPoint[0]
                lat = startPoint[1] + pointIndex*delta_lat
                #if Latitude is greater than 180 or smaller than -180, the loxodrome has reached the pole and must stop.
                if lat >= 180.0:
                    segment.append([lon,180.0])
                    return segment
                if lat <=-180.0:
                    segment.append([lon,-180.0])
                    return segment
                else:
                    segment.append([lon,lat])
        return segment

class BadArguments(Exception):
    def __init__(self, message):
        self.message=message
    def __str__(self):
        return self.message
