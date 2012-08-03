#!/usr/bin/env python

import glob, os, re, math, numpy, sys
import bmaudio
import pylab
import matplotlib
import numpy

from parameters import *
import headDistortion
import scipy.special 

sourceAzimuth = 90
f = spectralRange/float(spectrumBins) * numpy.arange(spectrumBins)
w = f*2*math.pi

print "R", R
print "r", r
print "c", c

degreesResolution=.05
degreesResolution=1
azimuths = numpy.arange(-90+.1,90,degreesResolution)
spectralBinPositions = numpy.arange(0,spectralRange, spectralRange/spectrumBins)
decodings=dict(
	inphase = headDistortion.decoding2dInPhase,
	maxre = headDistortion.decoding2dMaxRe,
	maxrv = headDistortion.decoding2dBasic,
	)
if len(sys.argv)<2 or sys.argv[1] not in decodings.keys() :
	print "Please specify one decoding:", ", ".join(decodings.keys())
	sys.exit()
if len(sys.argv)<3 or not int(sys.argv[2]) :
	print "Please specify an order"
	sys.exit()
maxOrder = int(sys.argv[2])

decodingName=sys.argv[1]
decoding=decodings[decodingName](maxOrder)
print "decoding", decoding, "order", maxOrder


print "Computing delays..."
distances = numpy.array([
	bmaudio.sphericalHeadDelay(numpy.radians(a),R,r)
	for a in azimuths
	])
relativeDistances = distances + R - r # force closer distance be 0
delays = relativeDistances / c  # distance -> time


print "Computing simplified delays..."
distancesSimplified = numpy.array([
	bmaudio.sphericalHeadDelaySimplified(numpy.radians(a),R,r)
	for a in azimuths
	])
relativeDistancesSimplified = distancesSimplified + R - r # force closer distance be 0
delaysSimplified = relativeDistancesSimplified / c  # distance -> time


print "Computing spherical harmonics..."
sphericalHarmonics = numpy.zeros((maxOrder+1, len(azimuths)))
for l in xrange(maxOrder+1) :
	sphericalHarmonics[l,:] = numpy.sqrt(2*l+1)/2. * scipy.special.legendre(l)(numpy.sin(numpy.radians(-azimuths)))


print "Computing components..."
components = numpy.zeros((maxOrder+1, spectrumBins), dtype=numpy.complex)
for l in xrange(maxOrder+1) :
	for a,sh,d in zip(azimuths, sphericalHarmonics[l,:], delaysSimplified) :
		components[l,:] += numpy.cos(numpy.radians(a))*sh*numpy.exp(-1j*w*d)

# TODO: review normalizations
components/=len(azimuths) # Normalize the discrete sum of len(azimuths) deltas
components*=numpy.pi
components/=components.max()
dbcomponents = 20*numpy.log10(numpy.abs(components))


print "Computing perceived filtering..."
filtering=numpy.zeros((len(azimuths), spectrumBins), dtype=numpy.complex)
for l in xrange(maxOrder+1) :
	for i,(a,sh) in enumerate(zip(azimuths,sphericalHarmonics[l,:])) :
		filtering[i,:] += decoding[l]*sh*components[l,:]
filteringdb = 20*numpy.log10(numpy.abs(filtering))




print "Computing components from time domain..."
componentsTD = numpy.zeros((maxOrder+1, nBins))
to=R/c
for l in xrange(maxOrder+1) :
	componentsTD[l,:] = numpy.sqrt(2*l+1)/2/to*numpy.array([
		scipy.special.legendre(l)(t/to) # -to,0
		if t<0 else (
		numpy.cos(t/to) * scipy.special.legendre(l)(numpy.sin(t/to))
		if t<to*numpy.pi/2 else
		0)
		for t in numpy.linspace(-to+.00000000001,nBins/float(samplingRate)-to,nBins)
		]) / nBins/numpy.pi/2
componentsFD = numpy.zeros((maxOrder+1, spectrumBins), dtype=numpy.complex)
for l in xrange(maxOrder+1) :
	componentsFD[l,:] = numpy.fft.rfft(componentsTD[l,:])/nBins/numpy.pi/2
componentsFD/=componentsFD.max()
dbcomponentsFD = 20*numpy.log10(numpy.abs(componentsFD))
print "Components FD max", componentsFD.max()


print "Computing percevied filtering from time domain analytic expression..."
filteringTD=numpy.zeros((len(azimuths), nBins))
filteringFD=numpy.zeros((len(azimuths), spectrumBins), dtype=numpy.complex)
for l in xrange(maxOrder+1) :
	for i,(a,sh) in enumerate(zip(azimuths,sphericalHarmonics[l,:])) :
		filteringTD[i,:] += decoding[l]*sh*componentsTD[l,:]
for i,(a,sh) in enumerate(zip(azimuths,sphericalHarmonics[l,:])) :
	filteringFD[i,:] = numpy.fft.rfft(filteringTD[i,:])
filteringFDdb = 20*numpy.log10(numpy.abs(filteringFD))

print filteringFDdb






def plotField(name, field, xaxis, yaxis) :
	levels = 30
	levels = xrange(-60,6,3)
	field[numpy.isinf(field)]=-6
	pylab.rcParams['figure.figsize']=(15,7)
	pylab.rcParams['figure.subplot.left'] = 0.08
	pylab.rcParams['figure.subplot.right'] = 1.0
	pylab.contourf(xaxis, yaxis, field, levels, cmap=pylab.cm.Greys_r, extend="both")
	pylab.colorbar().set_label("Magnitude (dB)")
	c=pylab.contour(xaxis, yaxis, field, levels, cmap=pylab.cm.Greys)
	pylab.ylim(-90,90)
	pylab.xlabel("Frequency (Hz)")
	pylab.ylabel("Source azimuth (degrees)")
	pylab.clabel(c)
	pylab.savefig('figures/'+os.path.basename(os.path.splitext(__file__)[0])+"-"+name+".pdf", format="pdf")
	pylab.savefig('figures/'+os.path.basename(os.path.splitext(__file__)[0])+"-"+name+".png", format="png")
#	pylab.show()
	pylab.close()

print "Plotting..."
if 0 :
	plotField(decodingName+"_%02d"%maxOrder, filteringdb, spectralBinPositions, azimuths)


else : # debug plots
#	plotField(decodingName, filteringTD[:,:50], [x for x in xrange(50)], azimuths)
#	pylab.plot(azimuths, delays, label="delay"); pylab.title("Delay in seconds of a sound comming from azimuth")
#	pylab.plot(azimuths, delaysSimplified, label="delaySimplified"); pylab.title("Delay in seconds of a sound comming from azimuth")
	for l in xrange(maxOrder+1) :
#		pylab.plot(azimuths, sphericalHarmonics[l], label=str(l)); pylab.title("Spherical Harmonics")
#		pylab.plot(f, numpy.abs(components[l,:]), label=str(l)); pylab.title("HRTF Components")
#		pylab.plot(numpy.abs(numpy.fft.irfft(components[l,:])), label=str(l)); pylab.title("HRTF Components (Time domain)")
#	for l in xrange(maxOrder+1) :
#		pylab.plot(componentsTD[l,:]+.2, label="td"+str(l)); pylab.title("HRTF Components in time domain"); pylab.ylim(-.5,.5); pylab.xlim(0,30)
		pylab.plot(f, dbcomponents[l,:], label="A"+str(l)); pylab.title("HRTF Components (Spectral Magnitude dB)"); pylab.ylim(-50,1)
#		pylab.plot(f, dbcomponentsFD[l,:], label="B"+str(l)); pylab.title("HRTF Components in time domain (Spectral Magnitude dB)"); pylab.ylim(-50,1)

	pylab.legend()
	pylab.show()



