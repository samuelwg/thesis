#!/usr/bin/env python

import glob, os, re, math, numpy, sys
import bmaudio
import pylab

from parameters import *

import headDistortion
sourceAzimuth = 90


decodings=dict(
	inphase = headDistortion.inphasePolarPattern,
	maxre = headDistortion.maxrePolarPattern,
	maxrv = headDistortion.basicPolarPattern,
	)
if len(sys.argv)<2 or sys.argv[1] not in decodings.keys() :
	print "Please specify one decoding:", ", ".join(decodings.keys())
	sys.exit()
decodingName=sys.argv[1]
useLog = '--log' in sys.argv
showPhase = '--phase' in sys.argv

polarPattern = decoding=decodings[decodingName]

plot = bmaudio.SpectrumDisplay()
plot.inDb()
if showPhase : plot.showPhase()
if useLog : plot.setLogFrequency()
plot.ylim(-25,5)
#plot.flim(0,20000)
#ordersToShow = [1,3,5,7,9,15,17]
#ordersToShow = numpy.arange(7)+32

cases = [
	("analytical", headDistortion.sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow)),
]
for legend, filters in cases :
	for order in ordersToShow :
		plot.addSpectrumData(
			filters[order], 22050, "Order %i"%(order))
plot.flim(0,spectralRange)
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+"-"+decodingName+".pdf", "pdf")




