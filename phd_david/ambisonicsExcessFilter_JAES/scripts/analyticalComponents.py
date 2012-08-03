#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion

def polarPatternX(azimuth, order) :
	return [ math.cos(order * azimuth) ]
def polarPatternY(azimuth, order) :
	return [ math.sin(order * azimuth) ]
N=721
evenOrdersToShow = 0,2,4#,6#,8,#22,24,34
oddOrdersToShow = 1,3,#5#,7,9,#23,25,35
filters = dict()
filters.update(headDistortion.sphericHeadSimulation(
	bmaudio.sphericalHeadDelay, 0, N, polarPatternX, evenOrdersToShow))
filters.update(headDistortion.sphericHeadSimulation(
	bmaudio.sphericalHeadDelay, 0, N, polarPatternY, oddOrdersToShow))
if True :
	plot = bmaudio.SpectrumDisplay()
	plot.inDb()
#	plot.showPhase()
	plot.ylim(-30,0)
	plot.flim(0,spectralRange)
	plot.setStylePreference(["lines","colors"])
#	plot.setStyleVariation("lines", ["-", "--"]) # even, odd
	plot.setStyleVariation("colors", ["k"])
	for order in sorted(filters.keys()) :
		print order,"\t",(abs(filters[order]).argmax())*spectralRange/nBins*R/c
	#	print order, filters[order].argmax()*spectralRange/nBins*R/c
		legend = "$\\mathrm{H_{%%s,%s}}$"% ('-' if order&1 else '+')
		plot.addSpectrumData( filters[order], spectralRange, legend%order)
#	plot.setLegendTitle("$\\mathrm{H_{l, \\sigma}}$ components")
	plot.legendPosition("upper right")
	

	plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")
else :
	x,y = zip(*[
		(order,(abs(filters[order]).argmax())*spectralRange/spectrumBins*R/c*2*numpy.pi)
		for order in sorted(filters.keys()) ])
	print y

	pylab.plot(x,y,"+-")
	pylab.plot([x+.5 for x in xrange(max(filters.keys()))],"--" )
	pylab.show()

"""
Contiguous order lobule distance is:
1=f_0*R*2*pi/c
f_0 = c/(2*pi*R)
w_0 = c/R
"""



