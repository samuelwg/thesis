#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *
import headDistortion

def figurePath(file, extension) :
	return os.path.join('figures', os.path.splitext(os.path.basename(file))[0] + "." + extension)

sourceAzimuth = 90

polarPattern = headDistortion.maxrePolarPattern
#polarPattern = headDistortion.maxrePolarPattern
ordersToShow = [0,1,2,3,5,7,9]
#ordersToShow = numpy.arange(7)+32
#ordersToShow = 1, 2, 5, 10
Ncontinuous = 720
continuousFilter = headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, Ncontinuous, polarPattern, ordersToShow)
maxN = 72
Ns = range(1,maxN)
progressions = dict([(order, numpy.zeros(maxN)) for order in ordersToShow ])
for N in Ns :
	filters = headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow)
	for order in ordersToShow :
		progressions[order][N] = 2* spectralRange / nBins * abs(bmaudio.spectrumMagnitudeInDb(filters[order]/continuousFilter[order])).cumsum().searchsorted(.4*nBins)

#for order, color in zip(ordersToShow, "bgrycmk"[:len(ordersToShow)]) :
#	pylab.plot(xrange(0,maxN,2), progressions[order][::2], "--"+color, label="O%i even"%order)
for order, color in zip(ordersToShow, ("bgrycmk"*4)[:len(ordersToShow)]) :
	pylab.plot(Ns[::2], progressions[order][::2], "-"+color, label="O%i odd"%order)
pylab.legend(loc=4)
pylab.show()

pylab.savefig(figurePath(__file__,"pdf"))



