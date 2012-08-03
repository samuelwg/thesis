#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab
import random


print "Simulating horizontal sources an spheric head..."

N = 40000
N = 72
r=1.4
R=0.088
nBins = 1024*4
samplingRate = 44100. # samples/s
spectralRange = samplingRate/2.
c = 344.
sphericalHarmonics = bmaudio.sphericalHarmonics2D

spectrumBins = nBins/2+1
shFilters = dict( [ (shName, numpy.zeros(spectrumBins, numpy.complex))
	for shOrder, shName, _, shFunction, shWeight in sphericalHarmonics ])

w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

print "Computing the filter for each spherical harmonic..."
for i in xrange(N) :
#	i = random.random() * N # Montecarlo sampling (instead of equidistant)
	e = 0.
	a = i * numpy.pi / N
	t = (bmaudio.sphericalHeadDelay(a+math.pi/2, R, r) - (r-R)) / c
	sinusoid = numpy.exp( 1j * t * w )
	for shOrder, shName, _, shFunction, shWeight in sphericalHarmonics :
		shFilters[shName] += shFunction(a,e) * sinusoid 

print "Obtaining the resulting filter for different orientations..."
steps = 72
orderFilters = {}
for azimuthDegrees in 360./steps*numpy.arange(steps) :
	orderFilters[azimuthDegrees] = {}
	for shOrder, shName, _, shFunction, shWeight in sphericalHarmonics :
		orderFilters[azimuthDegrees][shOrder] = numpy.zeros(spectrumBins, numpy.complex)
	for order in orderFilters[azimuthDegrees].keys() :
		for shOrder, shName, _, shFunction, shWeight in sphericalHarmonics :
			if shOrder>order : continue
			orderFilters[azimuthDegrees][order] += shWeight * shFunction(math.radians(azimuthDegrees),0) *  shFilters[shName]

#plot.showPhase()
#plot.ylim(0)
for azimuth, filters in orderFilters.iteritems() :
	plot = bmaudio.SpectrumDisplay()
	plot.inDb()
	for order, filter in filters.iteritems() :
		plot.addSpectrumData(filter, samplingRate/2, "Azimuth %s, Order %s"%(azimuth, order))
	plot.ylim(ymin=-10, ymax=50)
	plot.hardcopy(os.path.splitext(__file__)[0]+"%03i.png"%azimuth, "png")








