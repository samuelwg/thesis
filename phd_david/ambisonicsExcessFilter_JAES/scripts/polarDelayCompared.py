#!/usr/bin/env python
import math
import os
import bmaudio
import sys
import pylab
from parameters import *

def figurePath(file, extension) :
	return os.path.join('figures', os.path.splitext(os.path.basename(file))[0] + "." + extension)

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)

azimuths = numpy.arange(0,360,360/N)
radiansAzimuths = [math.radians(degrees) for degrees in azimuths]
sphericHeadDelays = dict([
	(azimuthDegrees, (
		bmaudio.sphericalHeadDelay(math.radians(azimuthDegrees), R, r) - (r-R)
		) / c
	)
	for azimuthDegrees in azimuths])
realHeadDelays = bmaudio.azimuthDelaysFromHorizontalHrtf(hrtfDatabase)
pylab.rcParams["figure.figsize"] = (8,5)
pylab.axes(polar=True).set_theta_zero_location("N")
pylab.axes(polar=True).set_thetagrids(
	range(0,360,45),
	labels= [str(g) for g in range(0,180+1,45)+range(-135,0,45)]
	)
pylab.polar(radiansAzimuths, numpy.array([sphericHeadDelays[azimuth] for azimuth in azimuths])*1000, "k--", label='Spherical head' )
pylab.polar(radiansAzimuths, numpy.array([realHeadDelays[azimuth] for azimuth in azimuths])*1000, "k-", label="Real head" )
pylab.ylabel("ms")
pylab.legend()
pylab.savefig(figurePath(__file__,"pdf"))
pylab.show()

