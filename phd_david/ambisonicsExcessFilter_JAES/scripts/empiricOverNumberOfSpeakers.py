#!/usr/bin/env python

import glob, os, re, math, numpy, sys
import bmaudio
import pylab
import headDistortion

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)

r=1.4
R=0.075 # interaural distance
R=0.088 # equivalent sphere
a0 = math.acos(R/r)
nBins = 1024*4
samplingRate = 44100
spectralRange = samplingRate/2.
c=344
spectrumBins = nBins/2+1
Ns = [72, 36, 24, 18, 12]

delays = bmaudio.azimuthDelaysFromHorizontalHrtf(hrtfDatabase)

w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)
filters = dict( [ (N, numpy.zeros(spectrumBins,numpy.complex)) for N in Ns ])
for N in Ns :
	for azimuth in xrange(0,360,360/N) :
		sinusoid = numpy.exp( 1j * delays[azimuth] * w )
		filters[N] += sinusoid /N * sum(headDistortion.maxrePolarPattern(math.radians(azimuth-90), 0))

plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.showPhase()
plot.ylim(-20,2)
plot.flim(fmax=20000)
for N in Ns :
	plot.addSpectrumData(filters[N], 20050, "%i speakers"%N)
plot.hardcopy(os.path.splitext(__file__)[0]+".pdf", "pdf")


