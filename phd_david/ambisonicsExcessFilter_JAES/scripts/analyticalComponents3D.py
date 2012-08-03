#!/usr/bin/env python


import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion

from parameters import *
samplingRate = 48000
nBins = 1024
spectrumBins = nBins/2+1

import scipy.special

ordersToShow = [0] + [ 1, 2, 3, 4, 5, 6, 7 ][:4]

def hrtfComponentsFromDeltas(ordersToShow, spectrumBins, spectralRange) :
	f = spectralRange/float(spectrumBins) * numpy.arange(spectrumBins)
	w = f*2*math.pi
	azimuths = numpy.arange(-90+.1,90,5)
	print "Computing simplified delays..."
	distancesSimplified = numpy.array([
		bmaudio.sphericalHeadDelaySimplified(numpy.radians(a),R,r)
		for a in azimuths
		])
	relativeDistancesSimplified = distancesSimplified + R - r # force closer distance be 0
	delaysSimplified = relativeDistancesSimplified / c  # distance -> time

	print "Computing 3D spherical harmonics..."
	sphericalHarmonics = numpy.zeros((len(ordersToShow), len(azimuths)))
	for i,l in enumerate(ordersToShow) :
		sphericalHarmonics[i,:] = numpy.sqrt(2*l+1)/2. * scipy.special.legendre(l)(numpy.sin(numpy.radians(-azimuths)))

	print "Computing components by adding freq domain deltas..."
	componentsFromDeltas = numpy.zeros((len(ordersToShow), spectrumBins), dtype=numpy.complex)
	for i,l in enumerate(ordersToShow) :
		for a,sh,d in zip(azimuths, sphericalHarmonics[i,:], delaysSimplified) :
			# cos(a) because is the sphere. the whole azimuth shares delay and gain
			componentsFromDeltas[i,:] += numpy.cos(numpy.radians(a))*sh*numpy.exp(-1j*w*d)
	componentsFromDeltas /= len(azimuths) # because we are adding those many deltas
	componentsFromDeltas *= numpy.pi # TODO: why? where we missed it

	return componentsFromDeltas


def hrtfComponentsFromTimeFourmula(ordersToShow, spectrumBins, spectralRange) :
	print "Computing time domain formula..."
	# Time domain formula
	componentTD = numpy.zeros((len(ordersToShow), nBins))
	for i,l in enumerate(ordersToShow) :
		componentTD[i,:] = headDistortion.sphericalHead3dSH(
			l, samplingRate, nBins)

	print "Obtaining freq domain by fft time domain..."
	componentsFromTD = numpy.zeros((len(ordersToShow), spectrumBins), dtype=numpy.complex)
	for i,l in enumerate(ordersToShow) :
		componentsFromTD[i,:] = numpy.fft.rfft(componentTD[i,:])
	componentsFromTD /= spectrumBins # because of the ifft requires so
	# TODO: Still missing a 95.4 factor!! so dumb normalization
	componentsFromTD /= componentsFromTD[0,0]

	return componentTD, componentsFromTD

def hrtfComponentsFromFreqFormula(ordersToShow, spectrumBins, spectralRange) :
	print "Computing freq domain formula..."
	# Freq domain formula
	componentsFromFD = numpy.zeros((len(ordersToShow), spectrumBins), dtype=numpy.complex)
	for i,l in enumerate(ordersToShow) :
		print l
		componentsFromFD[i,:] = headDistortion.sphericalHead3dSHSpectrum(l, spectralRange, spectrumBins)

	return componentsFromFD

componentsFromDeltas = componentsTD = componentsFromFD = componentsFromTD = None

componentsFromDeltas = hrtfComponentsFromDeltas(ordersToShow, spectrumBins, spectralRange)
componentsTD, componentsFromTD = hrtfComponentsFromTimeFourmula(ordersToShow, spectrumBins, spectralRange)
componentsFromFD = hrtfComponentsFromFreqFormula(ordersToShow, spectrumBins, spectralRange)

print componentsFromFD[0,0], componentsFromTD[0,0], componentsFromDeltas[0,0]

print "Plotting freq domain functions for the paper..."
plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.ylim(-30,0)
#	plot.showPhase()
legend = "$\\mathrm{H_{%s,0,+}}$"
plot.setStylePreference(["lines","colors"])
plot.setStyleVariation("colors", ["k"])
for i, order in enumerate(ordersToShow) :
	plot.addSpectrumData(componentsFromDeltas[i,:], spectralRange, legend%order)
plot.legendPosition("upper right")
plot.flim(0,spectralRange)
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")


if componentsTD is not None :
	print "Plotting time domain function..."
	for i,order in enumerate(ordersToShow) :
		pylab.plot(componentsTD[i,:], label=str(order))
	pylab.legend()
	pylab.show()

print "Plotting freq domain functions..."
plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.ylim(-30,5)
#	plot.showPhase()
legend = "$\\mathrm{H_{%s,0,+}}$"
if componentsFromTD is not None :
	for i, order in enumerate(ordersToShow) :
		plot.addSpectrumData(componentsFromTD[i,:], spectralRange, legend%order)
if componentsFromFD is not None :
	for i, order in enumerate(ordersToShow) :
		plot.addSpectrumData(componentsFromFD[i,:], spectralRange, legend%order)
if componentsFromDeltas is not None :
	for i, order in enumerate(ordersToShow) :
		plot.addSpectrumData(componentsFromDeltas[i,:], spectralRange, legend%order)
plot.legendPosition("upper right")
plot.flim(0,spectralRange)
plot.show()

print "Plotting components convolved by 0 order component inverse"
plot = bmaudio.SpectrumDisplay()
plot.inDb()
plot.ylim(-30,5)
#	plot.showPhase()
legend = "$\\mathrm{H_{%s,0,+}}$"
if componentsFromDeltas is not None :
	for i, order in enumerate(ordersToShow) :
		plot.addSpectrumData(componentsFromDeltas[i,:]/componentsFromDeltas[0,:], spectralRange, legend%order)
plot.legendPosition("upper right")
plot.flim(0,spectralRange)
plot.show()



