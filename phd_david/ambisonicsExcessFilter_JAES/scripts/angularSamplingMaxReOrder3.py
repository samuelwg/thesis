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
plot.ylim(-20,15)
plot.setStylePreference(["lines","colors"])
#plot.setStyleVariation("colors",["k"])
#plot.flim(0,20000)
order = 200
ordersToShow = order,
Ncontinuous = 722
continuousFilter = headDistortion.sphericHeadSimulation(
	bmaudio.sphericalHeadDelay,
	sourceAzimuth,
	Ncontinuous,
	polarPattern,
	ordersToShow) [order]
Ns = 6, 12, 18, 24, 30, 36, 42, 45, 48, 72,
Ns = xrange(3,70,2)

filters = [
	headDistortion.sphericHeadSimulation(
		bmaudio.sphericalHeadDelay,
		sourceAzimuth,
		N,
		polarPattern,
		ordersToShow)
		[order]
	for N in Ns ]

divergences = [
	filter/continuousFilter
	for filter in filters ]

divergencePoints = [
	(numpy.abs(numpy.abs(divergence)-1)<.02).sum()
	for divergence in divergences ]

pylab.plot(Ns,divergencePoints,'x-')
pylab.grid()
pylab.show()

for N, divergence in zip(Ns,divergences) :
	plot.addSpectrumData(
		divergence, 22050, "%i speakers"%(N))

plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")

"""
Conclusions:
Odd number of speakers converge faster than even numbers
Divergence point is linear with the number of speakers (considering just event or just odd)
The velocity of convergence is faster as we rise the
The decoding criterion has influence on the 
maxre nspeakers<=order do not match if nspeakers even
maxre nspeakers<order/2 do not match if nspeakers odd

"""



