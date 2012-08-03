#!/usr/bin/env python

import glob, os, re, math, numpy, sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
from bmaudio import *
import pylab

databaseFile = selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = HrtfDatabase(databaseFile)

deviations=[]
sampleDelaysPI=[]
sampleDelaysMP=[]
sampleDelays85=[]
azimutalPI=[]
azimutal85=[]
azimutalMP=[]
degrees=[]

for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
	from math import radians
	elevation = radians(elevationDegrees)
	azimuth = radians(azimuthDegrees)
	samplingRate, data = loadWave(response)
	sampleDelaysPI.append(delayWithMinimumPhaseCorrelation(data,filterFreq=20000, interpolate=True))
	if elevation == 0 :
		azimutalPI.append(sampleDelaysPI[-1])
		degrees.append(azimuthDegrees)


r=1.4
R=0.075 # interaural distance
R=0.088 # equivalent sphere
a0 = math.acos(R/r)
c=344

if True :
	radians = numpy.array([math.radians(a) for a in degrees])
	analytic = numpy.array([samplingRate/c*(sphericalHeadDelay(math.radians(a), R, r)) for a in degrees])
	pylab.polar(radians, numpy.array(azimutalPI)-min(azimutalPI), label='HRTF delay (empirical)')
	pylab.polar(radians, analytic-min(analytic), label='HRTF delay (analytic)')
	pylab.title("Analytic spherical head vs. KEMAR (horizontal plane)")
	pylab.legend()
	pylab.show()

if True :
	pylab.title("HRTF Delay histogram")
	for label, sampleDelays in (
		("Number of points",sampleDelaysPI),
		) :
		histogram,edges=numpy.histogram(sampleDelays,range=(0,512),bins=512)
		pylab.plot(histogram[0:100], label=label)
	pylab.legend()
	pylab.show()

minPhase = minimumPhase(histogram)
fft=numpy.fft.rfft(histogram)
fftMinPhase = numpy.fft.rfft(minPhase)

if True :
	plot = SpectrumDisplay()
	plot.inDb()
	plot.showPhase()
	plot.addSpectrumData(fft,samplingRate/2, "Original")
	plot.addSpectrumData(fftMinPhase, samplingRate/2, "Min Phase")
	plot.addSpectrumData(1/fftMinPhase, samplingRate/2, "Inverse")
	plot.show()

headCompensation = numpy.fft.irfft(fft)
headCompensationMinPhase = numpy.fft.irfft(fftMinPhase)
headCompensationMinPhaseFilter = numpy.fft.irfft(1/removeZeros(fftMinPhase))

if False :
	pylab.title("Time domain filters")
	pylab.plot(headCompensation, label="Original")
	pylab.plot(headCompensationMinPhase, label="Min Phase")
	pylab.plot(headCompensationMinPhaseFilter, label="Inverse")
	pylab.legend()
	pylab.show()


def normalize(data) :
	factor = max(abs(data))
	data /= factor
	return factor


print "headCompensation norm factor: ", normalize(headCompensation)
print "headCompensationMinPhase norm factor: ", normalize(headCompensationMinPhase)
print "headCompensationMinPhaseFilter norm factor: ", normalize(headCompensationMinPhaseFilter)

if False :
	pylab.title("Time domain filters")
	pylab.plot(headCompensation, label="Original")
	pylab.plot(headCompensationMinPhase, label="Min Phase")
	pylab.plot(headCompensationMinPhaseFilter, label="Inverse")
	pylab.legend()
	pylab.show()

saveWave("headCompensationMinPhase.wav", headCompensationMinPhase, samplingRate)
saveWave("headCompensation.wav", headCompensation, samplingRate)
saveWave("headCompensationFilter.wav", headCompensationMinPhaseFilter, samplingRate)





