#!/usr/bin/env python
"""\
Creates a polar plot with the gains of a 2d ambisonics decoding of a plane wave
at a given order and decoding criteria.

Usage: ./flatWaveDecoding2d.py maxre|maxrv|inphase [order1 [ order2 [ order3 ... ] ] ]
"""

import math
import os
import bmaudio
import pylab
import sys
from parameters import *
import headDistortion

decodings=dict(
	inphase = headDistortion.inphasePolarPattern,
	maxre = headDistortion.maxrePolarPattern,
	maxrv = headDistortion.basicPolarPattern,
	)
if len(sys.argv)<2 or sys.argv[1] not in decodings.keys() :
	print "Please specify one decoding:", ", ".join(decodings.keys())
	sys.exit()
decodingName=sys.argv[1]
polarPattern = decoding=decodings[decodingName]
title = dict(
	maxre = "Max $r_E$",
	maxrv = "Max $r_V$",
	inphase = "In-phase",
	)[decodingName]

ordersToShow = [1,3,9,35]
if  sys.argv[-1:1:-1] : ordersToShow = [int(a) for a in sys.argv[2:]]

def figurePath(file, infix, extension) :
	return os.path.join('figures', os.path.splitext(os.path.basename(file))[0] + infix + "." + extension)

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
pylab.title(title,horizontalalignment='center', verticalalignment='baseline', position=(.5,-.13))
pylab.savefig(figurePath(__file__,"-"+decodingName,"pdf"))
pylab.show()


