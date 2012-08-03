#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion
sourceAzimuth = 90

polarPattern = headDistortion.basicPolarPattern
plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.showPhase()
plot.ylim(-25,5)
#plot.flim(0,20000)
#ordersToShow = [1,3,5,7,9,15,17]
#ordersToShow = numpy.arange(7)+32

cases = [
	("measured", headDistortion.hrtfFullSimulation(headDistortion.hrtfDatabase, sourceAzimuth, N, polarPattern, ordersToShow)),
	("delay", headDistortion.hrtfDelaySimulation(headDistortion.hrtfDatabase, sourceAzimuth, N, polarPattern, ordersToShow)),
	("analytical", headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow)),
]
for legend, filters in cases :
	for order in ordersToShow :
		plot.addSpectrumData(
			filters[order], 22050, "%s o%i"%(legend,order))
plot.flim(0,21000)
plot.setStyleLegend(
	"Order", ordersToShow,
	"Approximation", [legend for legend, filters in cases])
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")




