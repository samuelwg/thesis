#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion
sourceAzimuth = 90

polarPattern = headDistortion.inphasePolarPattern
plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.showPhase()
plot.ylim(-25,5)
#plot.flim(0,20000)
#ordersToShow = [1,3,5,7,9,15,17]
#ordersToShow = numpy.arange(7)+32
for legend, filters in dict(
	RealDelay = headDistortion.hrtfDelaySimulation(headDistortion.hrtfDatabase, sourceAzimuth, N, polarPattern, ordersToShow),
	Spheric = headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow),
).items() :
	for order in ordersToShow :
		plot.addSpectrumData(
			filters[order], 22050, "%s o%i"%(legend,order))
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")




