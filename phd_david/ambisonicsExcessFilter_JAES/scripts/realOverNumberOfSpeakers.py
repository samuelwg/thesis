#!/usr/bin/env python

import glob, os, re, math, numpy, sys
import bmaudio
import pylab
from parameters import *

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)

Ns = [72, 36, 24, 18, 12]
nBins = 1024*4

spectrumBins = nBins/2+1
filters = dict( [ (N, numpy.zeros(spectrumBins,numpy.complex)) for N in Ns ])
w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

hrtfs = {}
delays = {}

for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
	if elevationDegrees != 0 : continue
	from math import radians
	e = radians(elevationDegrees)
	a = radians(azimuthDegrees)
	samplingRate, data = bmaudio.loadWave(response)
	t = bmaudio.delayWithMinimumPhaseCorrelation(data,filterFreq=20000, interpolate=True) / samplingRate
	delays[azimuthDegrees] = t
	hrtfs[azimuthDegrees] = numpy.fft.rfft(data,nBins)


minDelay = min(delays.values())

for N in Ns :
	for azimuth in xrange(0,360,360/N) :
		filters[N] += hrtfs[azimuth] / N
for N in Ns :
	filters[N] /= hrtfs[90]

plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.showPhase()
plot.ylim(-40)
for N in Ns :
	plot.addSpectrumData(filters[N][0:1500], samplingRate/2, "%i speakers"%N)
plot.hardcopy(os.path.splitext(__file__)[0]+".pdf", "pdf")


