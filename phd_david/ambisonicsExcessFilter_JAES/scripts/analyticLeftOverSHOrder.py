#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab
import random


print "This scripts simulate the resulting on a spherical head"
print "given several layouts of speakers and an infinite layout (Montecarlo aproximation)."

from parameters import *

analitic = dict( [ (order, numpy.zeros(spectrumBins,numpy.complex)) for order in ordersToShow])

w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

print "Simulating delays from a left  in the horizontal plane..."
for i in xrange(N) :
	e = 0.
	a = i * 2* numpy.pi / N
	t = (bmaudio.sphericalHeadDelay(a+math.pi/2, R, r) - (r-R)) / c
	sinusoid = numpy.exp( 1j * t * w )
	ncos = numpy.array( [math.cos(n*a) for n in xrange(NOrder)] )
	for order in ordersToShow :
		analitic[order] += sum(ncos[0:order+1]) * sinusoid

for key in analitic.iterkeys() :
	analitic[key] /= numpy.max(analitic[key])

plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.showPhase()
plot.ylim(-25,2)
for order in ordersToShow :
	plot.addSpectrumData(analitic[order][:int(20000.*spectrumBins*2/samplingRate)], 20000, "Order %i"%order)
plot.hardcopy(os.path.splitext(__file__)[0]+".pdf", "pdf")



