#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab
import random
from parameters import  *

print "This scripts simulate the computation of the W component on a spherical head"
print "given several layouts of speakers and an infinite layout (Montecarlo aproximation)."

Ns = [72, 36, 24, 18, 12]
nBins = 1024*4
samplingRate = 44100. # samples/s
spectralRange = samplingRate/2.
c = 344.
sphericalHarmonics = bmaudio.sphericalHarmonics[0:1] # just take the W component

spectrumBins = nBins/2+1
filters = dict( [ (str(N), numpy.zeros(spectrumBins,numpy.complex)) for N in Ns ])

w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

print "Simulating delays from %s equiangular speakers in the horizontal plane..."%(", ".join(str(i) for i in Ns))
for N in Ns : 
	for i in xrange(N) :
		e = 0.
		a = i * numpy.pi / N
		t = (bmaudio.sphericalHeadDelay(a+math.pi/2, R, r) - (r-R)) / c
		sinusoid = numpy.exp( 1j * t * w )
		shOrder, shName, _, shFunction, shWeight = sphericalHarmonics[0]
		filters[str(N)] += (shWeight * shFunction(a,e) / N) * sinusoid 
"""
N = 40000
print "Simulating delays from %i stochastic points in the horizontal plane..."%N
filters['stochastic'] = numpy.zeros(spectrumBins, numpy.complex)
for i in xrange(N) :
	i = random.random() * N # Montecarlo sampling (instead of equidistant)
	e = 0.
	a = i * numpy.pi / N
	t = (bmaudio.sphericalHeadDelay(a+math.pi/2, R, r) - (r-R)) * samplingRate / c
	sinusoid = numpy.exp( 1j * t * w )
	shOrder, shName, _, shFunction, shWeight = sphericalHarmonics[0]
	filters['stochastic'] += shWeight * shFunction(a,e) / N * sinusoid 
"""

sortedKeys = filters.keys()
sortedKeys.sort()

plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.showPhase()
plot.ylim(-40)
for name in sortedKeys :
	plot.addSpectrumData(filters[name], samplingRate/2, name)
plot.hardcopy(os.path.splitext(__file__)[0]+".pdf", "pdf")



