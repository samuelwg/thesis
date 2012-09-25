#!/usr/bin/python

import sys
#import scipy.io
import os
#import pylab
import numpy
import math
import wave

normalizationFactor = 4.
delay = 60 # in samples

def saveWave(filename, data, samplerate, verbose=False) :
	import wavefile
	with wavefile.WaveWriter(
		filename,
		channels=1,
		samplerate=samplerate,
		format = wavefile.Format.WAV | wavefile.Format.FLOAT
		) as writer :
		writer.write(data[:,numpy.newaxis])

def loadWave(filename, verbose=False) :
	import wavefile
	with wavefile.WaveReader(filename) as reader :
		data = numpy.zeros((reader.frames, reader.channels))
		return reader.samplerate, data

outputPattern = "kreuzer2/output_e%+03.f_a%05.1f.wav"
print "Reading positions..."
directions = numpy.loadtxt("octave/Kreuzer_Database/pairs.mat")
print "Reading magnitudes..."
impulseResponses = numpy.loadtxt("octave/Kreuzer_Database/hrtfsample.mat").transpose()
print "Reading phases..."
phases = numpy.loadtxt("octave/Kreuzer_Database/Phase.mat").transpose()

print "Generating files..."

maxError = 0
maxFirstBinError = 0
os.system('mkdir -p kreuzer2')
hrtfsFile = open('kreuzerDatabase2.hrtfs','w')
for (elevation, azimuth), magnitude, phase in zip(directions, impulseResponses, phases) :
	elevation = int(elevation)
	azimuth = int(azimuth)
	wavefile = outputPattern%(elevation,azimuth)
	print >> hrtfsFile, elevation, azimuth, wavefile
	spectrum = 10**(magnitude/20)*numpy.exp(-1j*math.pi/180*phase) # the phase is given inverted!
	# Zero padding to get 48KHz instead of 40KHz
	spectrum.resize(241)
	# Adding the 0 Hz bin
	spectrum = numpy.roll(spectrum,1)
	spectrum[0]=1
	# Adding a delay of 'delay' samples (TODO: copied from previous db code, not tested!)
	nbins=len(spectrum)
	k1=-delay*math.pi/nbins
	spectrum *= numpy.exp(k1*1j*numpy.arange(0,nbins))

	audio = numpy.fft.irfft(spectrum) / normalizationFactor
	saveWave(wavefile, audio, 48000)
	samplerate2, savedAudio = loadWave(wavefile)
	loadedSpectrum = numpy.fft.rfft(savedAudio)
	maxFirstBinError = max(maxFirstBinError,abs(normalizationFactor*loadedSpectrum[0]-spectrum[0])/abs(spectrum[0]))
	maxError = max(maxError,math.sqrt(numpy.mean((audio-savedAudio)**2)))
	if max(abs(audio))>1. :
		print >> sys.stderr, "Normalization factor of %s still produces clipping in %s. Needed %f"%(
			normalizationFactor, wavefile, normalizationFactor*max(abs(audio)))
		sys.exit(-1)
#	pylab.plot(audio)
#	pylab.plot(numpy.arctan2(spectrum.real, spectrum.imag))
	phase = numpy.arctan2(spectrum.real, spectrum.imag)
#	if azimuth == 40 :
#		pylab.plot(abs(20*numpy.log10(abs(spectrum[:201]))))
#		print (abs(20*numpy.log10(abs(spectrum[:201]))))
#		pylab.plot((180/math.pi*(phase[1:]-phase[:-1]+math.pi)-50)%360.-180.)
hrtfsFile.close()
#pylab.show()
print "The database was prepared OK, because the quantization error was only ", maxError, maxFirstBinError

