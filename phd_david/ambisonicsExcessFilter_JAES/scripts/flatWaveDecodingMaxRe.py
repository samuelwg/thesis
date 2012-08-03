#!/usr/bin/env python
import math
import os
import bmaudio
import pylab
from parameters import *
import headDistortion

ordersToShow=[1,3,9,35]

def figurePath(file, extension) :
	return os.path.join('figures', os.path.splitext(os.path.basename(file))[0] + "." + extension)

polarPattern = headDistortion.maxrePolarPattern
N*=2
azimuths = numpy.arange(0,361,360/N)
radiansAzimuths = [math.radians(degrees) for degrees in azimuths]
pylab.rcParams['figure.figsize']=(5,5)
for order in ordersToShow :
	paternValues = numpy.array([sum(polarPattern(azimuth, order))/(order+1) for azimuth in radiansAzimuths])
	paternValues /= max(paternValues)
	pylab.polar(radiansAzimuths, paternValues, label='Order %s'%order )
pylab.rgrids(numpy.arange(.4,1,.2),angle=220)
pylab.legend(loc=2)
pylab.title("Max Re",horizontalalignment='center', verticalalignment='baseline', position=(.5,-.13))
pylab.savefig(figurePath(__file__,"pdf"))
pylab.show()


