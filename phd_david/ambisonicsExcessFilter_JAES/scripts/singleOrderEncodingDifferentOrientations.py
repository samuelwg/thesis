#!/usr/bin/env python

import glob, os, re, math, numpy, sys
import bmaudio
import pylab
import matplotlib

from parameters import *
import headDistortion
sourceAzimuth = 90
maxOrder = 3
polarPattern = headDistortion.basicPolarPattern

def plotAllAngles(name, polarPattern, order, degreesResolution=5) :
		levels=30
		azimuths = xrange(-90,90,degreesResolution)
		spectralBinPositions = numpy.arange(0,spectralRange, spectralRange/spectrumBins)
		spectrums = numpy.zeros((len(azimuths), spectrumBins))
		pylab.rcParams['figure.figsize']=(15,7)
		pylab.rcParams['figure.subplot.left'] = 0.08
		pylab.rcParams['figure.subplot.right'] = 1.0
		for i,azimuth in enumerate(azimuths) :
			spectrums[i,:] = 20*numpy.log(numpy.abs(headDistortion.sphericHeadSimulation(
					bmaudio.sphericalHeadDelay, azimuth, N, polarPattern, [order]).values()[0]))
		pylab.contourf(spectralBinPositions, azimuths, spectrums,levels,cmap=pylab.cm.Greys_r)
		pylab.colorbar().set_label("Magnitude (dB)")
		pylab.contour(spectralBinPositions, azimuths, spectrums,levels,cmap=pylab.cm.Greys)
		pylab.xlabel("Frequency (Hz)")
		pylab.ylabel("Azimuth (degrees)")
		pylab.savefig('figures/'+os.path.basename(os.path.splitext(__file__)[0])+"-"+name+".pdf", format="pdf")
		pylab.show()
		pylab.close()

for name, decoding in [
	("inphase", headDistortion.inphasePolarPattern),
#	("maxrv", headDistortion.basicPolarPattern),
#	("maxre", headDistortion.maxrePolarPattern),
	] :
	plotAllAngles(name, decoding, order=maxOrder)



sys.exit()


plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.showPhase()
plot.ylim(-25,5)
#plot.flim(0,20000)
#ordersToShow = [1,3,5,7,9,15,17]
#ordersToShow = numpy.arange(7)+32
plot.setStylePreference(["lines","colors"])
plot.setStyleVariation("colors",["k"])
for legend, filters in [
	(label, 
	headDistortion.sphericHeadSimulation(
		bmaudio.sphericalHeadDelay, sourceAzimuth, N, polarPattern, ordersToShow))
	for label, sourceAzimuth in [
		("Ipsilateral", 90),
		("Contralateral", -90),
		("Frontal/Back", 0),
		]
	] :
	for order in ordersToShow :
		plot.addSpectrumData(
			filters[order], 22050, "%s"%legend)
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")




