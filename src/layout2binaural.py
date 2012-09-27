#!/usr/bin/python

# This scripts requires ipyclam python module installed

import ipyclam
import audio3d.hrtf

n=ipyclam.Network(ipyclam.Clam_NetworkProxy())

def addChannel(n, channel, leftIR, rightIR) :
	channelName = "_%02d"%channel

	AudioWindowing = "AudioWindowing"+channelName
	Convolution_L = "Convolution_L"+channelName
	Convolution_R = "Convolution_R"+channelName
	HrtfLoader_L = "HrtfLoader_L"+channelName
	HrtfLoader_R = "HrtfLoader_R"+channelName
	IFFT_L = "IFFT_L"+channelName
	IFFT_R = "IFFT_R"+channelName
	L = "L"+channelName
	R = "R"+channelName
	MyFFT = 'MyFFT'+channelName
	toStream_L = 'toStream_L'+channelName
	toStream_R = 'toStream_R'+channelName

	n [AudioWindowing] = 'AudioStream2Buffer'
	n [Convolution_L] = 'LowLatencyConvolution'
	n [Convolution_R] = 'LowLatencyConvolution'
	n [HrtfLoader_L] = 'ImpulseResponseLoader'
	n [HrtfLoader_L].ImpulseResponse = leftIR
	n [HrtfLoader_R] = 'ImpulseResponseLoader'
	n [HrtfLoader_R].ImpulseResponse = rightIR

	n [IFFT_L] = 'MyIFFT'
	n [IFFT_R] = 'MyIFFT'
	n [MyFFT] = 'MyFFT'
	n [toStream_L] = 'AudioBuffer2Stream'
	n [toStream_R] = 'AudioBuffer2Stream'

	n.Inputs[str(channel)] > n[AudioWindowing]["Audio stream"]
	n[AudioWindowing]["Audio buffer"] > n[MyFFT]["Audio Buffer"]
	n[Convolution_L].Output > n[IFFT_L]["Complex Spectrum"]
	n[Convolution_R].Output > n[IFFT_R]["Complex Spectrum"]
	n[HrtfLoader_L].W > n[Convolution_L].ImpulseResponse
	n[HrtfLoader_R].W > n[Convolution_R].ImpulseResponse
	n[IFFT_L]["Audio Buffer"] > n[toStream_L]["Audio buffer"]
	n[IFFT_R]["Audio Buffer"] > n[toStream_R]["Audio buffer"]
	n[MyFFT]["Complex Spectrum"] > n[Convolution_L].Input
	n[MyFFT]["Complex Spectrum"] > n[Convolution_R].Input
	n[toStream_L]["Audio stream"] > n.Adder_L["Input %i"%(channel-1)]
	n[toStream_R]["Audio stream"] > n.Adder_R["Input %i"%(channel-1)]

def bounce2binauralNetwork(leftIR, rightIR) :
	n.Inputs = "AudioSource"
	n.Inputs.NSources = len(leftIR) + 1 # for the subwoofer
	n.Adder_L = "AudioMixer"
	n.Adder_L.NumberOfInPorts = len(leftIR)
	n.Adder_R = "AudioMixer"
	n.Adder_R.NumberOfInPorts = len(leftIR)
	n.L = 'AudioSink'
	n.R = 'AudioSink'
	n.Adder_L > n.L
	n.Adder_R > n.R

	for i, (left, right) in enumerate(zip(leftIR, rightIR)) :
		addChannel(n, i+1, left, right)


import bmaudio
import sys

layoutfile = sys.argv[1]

database = audio3d.hrtf.selectHrtfDatabase(sys.argv)
hrtf = audio3d.hrtf.HrtfDatabase(database)

layout = []
for line in file(layoutfile) :
	if "#" in line : line = line[:line.index('#')]
	if not line.strip() : continue
	e,a,label = line.split()
	layout.append(( float(e), float(a)))


leftIR = [
	hrtf.nearestWavefile(e,+a)
	for e,a in layout
]
rightIR = [
	hrtf.nearestWavefile(e,-a)
	for e,a in layout
]

bounce2binauralNetwork(leftIR, rightIR)

sys.stdout.write(n.xml())






