#!/usr/bin/env python

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
import bmaudio
import Gnuplot
import numpy
import math
import re

def fileToSpectrum(file) :
	samplerate, wave = bmaudio.loadWave(file)
	return numpy.fft.rfft(wave)

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)
Iw = [ fileToSpectrum(file) for elevationDegrees, azimuthDegrees, file in hrtfDatabase._data ]
Iwsum = sum(abs(iw) for iw in Iw)/len(Iw)

Ew = fileToSpectrum(hrtfDatabase.equivalentPath("w"))
Ex = fileToSpectrum(hrtfDatabase.equivalentPath("x"))
Ey = fileToSpectrum(hrtfDatabase.equivalentPath("y"))
Ez = fileToSpectrum(hrtfDatabase.equivalentPath("z"))
nBins = len(Ew)
c=343
standardDistanceError=0.0067
# exp( (sigma*w)^2/2); sigma=d/c; w=i*spectralRange*pi/nBins/180
a=standardDistanceError*2*math.pi*22050/nBins/c
exponential = numpy.array([math.exp(-(a*n)**2/2) for n in xrange(0,nBins)])
normFactor = sum(1/exponential)

gp = bmaudio.SpectrumDisplay()

filter = abs(numpy.fft.irfft(1/exponential))[:nBins]
gp.inDb()
for spectrum, name in [
	(numpy.zeros(nBins), ""),
	(abs(Ew), "Ew"),
	(abs(Ex), "Ex"),
	(abs(Ey), "Ey"),
	(abs(Ez), "Ez"),
	(abs(Ew)/exponential, "Fixed Ew"),
	(abs(Ex)/exponential, "Fixed Ex"),
	(abs(Ey)/exponential, "Fixed Ey"),
	(abs(Ez)/exponential, "Fixed Ez"),
	(Iwsum, "Iwsum"),
	(abs(Ew)/Iwsum, "Ew/Iwsum"), # The ratio
#	[bin.real for bin in spectrum],
#	[bin.imag for bin in spectrum],
	(exponential, "Fix filter"),
	] :
	gp.addSpectrumData(spectrum, 256, name)
gp.show()
gp.hardcopy("IncoherenceSimulation.png","png") 
print filter/max(filter)
bmaudio.saveWave("coherenceFilter.wav", filter/max(filter), 44100)


#bmaudio.saveWave("lala.wav",computedInverseAudio,44100)

