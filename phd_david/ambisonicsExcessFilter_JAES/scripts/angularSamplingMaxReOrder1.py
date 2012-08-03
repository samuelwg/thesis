#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion
sourceAzimuth = 90

polarPattern = headDistortion.maxrePolarPattern
plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.showPhase()
plot.ylim(-15,15)
plot.flim(0,spectralRange)
plot.setStylePreference(["lines","colors"])
plot.setStyleVariation("colors",["k"])
plot.setLogFrequency()
plot.legendPosition("lower left")
#plot.flim(0,20000)
#ordersToShow = [1,3,5,7,9,15,17]
#ordersToShow = numpy.arange(7)+32
ordersToShow = 1,
Ncontinuous = 720
continuousFilter = headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, Ncontinuous, polarPattern, ordersToShow)
Ns = 12, 24, 36, 72,
for N in Ns :
	filters = headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow)
	for order in ordersToShow :
		print numpy.sum(numpy.abs(filters[order])**2), numpy.sum(numpy.abs(filters[order]**2))*math.sqrt(N)
		plot.addSpectrumData(
			filters[order]/continuousFilter[order], 22050, "%i speakers"%(N))
#			filters[order], 22050, "%i speakers"%(N))

plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")



