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
	wavefile = wave.open( filename, 'w' )
	wavefile.setframerate( samplerate )
	wavefile.setnchannels( 1 )
	wavefile.setsampwidth( 2 ) # number of bytes per sample
	wavefile.writeframesraw((32767.0*data).astype(numpy.int16).tostring())
	wavefile.close()

def loadWave(filename, verbose=False) :
	wavefile = wave.open( filename, 'r' )
	samplingRate   = wavefile.getframerate()
	bytesPerSample = wavefile.getsampwidth()
	nFrames        = wavefile.getnframes()
	data           = wavefile.readframes( nFrames )
	wavefile.close()
	if bytesPerSample==2:  # 16 bits per sample
		data = numpy.fromstring( data, numpy.int16 ) / 32767.0 # -1..1 values
	elif bytesPerSample==1: # 8 bits per sample
		data = (numpy.fromstring( data, numpy.uint8 ) / 128.0 ) - 1.0 # -1..1 values
	elif bytesPerSample==4 : # 32 bits per sample
		data = (numpy.fromstring( data, numpy.int32 ) / 2.**(32-1) ) # -1..1 values
	else :
		print bytesPerSample
		2+"asdfasdf"
	return samplingRate, data

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

