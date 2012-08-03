#!/usr/bin/env python

r=5 # mean distance from the speakers to the microphone in meters
standardMicrophoneDeviation = 0.005 # position deviation in times r
nSpeakers = 500
nBins = 512
sampleRate = 44100.

c = 300. # meters per second
spectralRange = sampleRate / 2

import numpy
import math
import RandomArray
import stats


T=r/c


spectrum = numpy.zeros( nBins-1, numpy.complex)
for x in RandomArray.normal(0,standardMicrophoneDeviation,(nSpeakers,)) :
	t1root = 1+x
	t2root = 1-x
	print t1root, t2root
	spectrum += numpy.array([
		complex(1, t1root*w*T) * math.e**complex(math.cos(w*T*t1root), math.sin(w*T*t1root)) /x / w / w / T / T -
		complex(1, t2root*w*T) * math.e**complex(math.cos(w*T*t2root), math.sin(w*T*t2root)) /x / w / w / T / T 
		for w in [ math.pi*2*normfreq*spectralRange/nBins 
			for normfreq in xrange(1,nBins) ] ])


import Gnuplot
gp=Gnuplot.Gnuplot(persist=1)
gp('set data style lines')
gp.plot(abs(spectrum), [bin.real for bin in spectrum], [bin.imag for bin in spectrum], numpy.zeros(nBins))
gp.hardcopy(filename="IncoherenceSimulation.png",terminal="png") 



speakerDistances = RandomArray.normal(0,standardMicrophoneDeviation,(nSpeakers,))
print speakerDistances
#print stats.stdev(speakerDistances), standardMicrophoneDeviation, stats.mean(speakerDistances)

spectrum = numpy.zeros( nBins, numpy.complex)
for x in speakerDistances :
	for gamma in range(0,360,1) :
		gammaRadiants = gamma*math.pi/180.
		sinGamma = math.sin(gammaRadiants)
		cosGamma = math.cos(gammaRadiants)
		for normfreq in xrange(nBins) :	
			w = math.pi*2*normfreq*spectralRange/nBins
			displacement = r*math.sqrt(1+2*x*cosGamma+x**2)  -r
			phase = w*displacement/c
			spectrum[normfreq] += sinGamma * complex(math.cos(phase), math.sin(phase) )
inverse = nSpeakers/spectrum

import Gnuplot
gp=Gnuplot.Gnuplot(persist=1)
gp('set data style lines')
gp.plot(abs(spectrum), [bin.real for bin in spectrum], [bin.imag for bin in spectrum], numpy.zeros(nBins))
gp.hardcopy(filename="IncoherenceSimulation.png",terminal="png") 




"""
import sympy
import sympy.numerics

t=sympy.Symbol("t")
mu=sympy.Symbol("mu")
sigma=sympy.Symbol("s")

normal=1./(sigma*math.sqrt(2*math.pi))/math.e**((t**2)/(2*sigma**2))
print "normal", normal

sine=math.e**(sympy.numerics.ComplexFloat(0, 2*math.pi*mu*t/nBins))
print "sine", sine

print "product", normal*sine

integral = sympy.integrate(normal*sine,t)
print "integral", integral

sympy.Plot(integral.subs(sigma,standardDelay*1))
"""


