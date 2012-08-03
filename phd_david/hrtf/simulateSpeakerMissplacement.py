#!/usr/bin/env python

print """\
This script compares the simulated impulse response of 
several incoherent sources following a normal distribution
on the delay and the theoretical response derived from the
analytical expression given an infinite number of speakers.
It also displays their corresponding compensating filters.
Adjustable parameters are on the beginning of the script.
"""

standardSpeakerDeviation = .016 # position deviation in meters
nSpeakers = 828
nBins = 256
sampleRate = 44100.
profile = False

c = 343 # meters per second


spectralRange = sampleRate / 2
standardDelay = standardSpeakerDeviation / c

def main() :

	import numpy
	import math
	import sys, os
	import RandomArray
	sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
	import bmaudio

	# Simulation of random error
	speakerDistances = numpy.random.normal(loc=0,scale=standardDelay,size=nSpeakers)

	simulatedSpectrum = numpy.array( [
		sum( complex(math.cos(w*delay), math.sin(w*delay)) for delay in speakerDistances )
			for w in ( math.pi*2*wavenumber*spectralRange/nBins for wavenumber in xrange(nBins) )
		])
	simulatedSpectrum/=nSpeakers
	simulatedInverse = 1/simulatedSpectrum
	simulatedInverseAudio = numpy.fft.irfft(simulatedInverse)

	# Analytical computation
	computedSpectrum = numpy.array([
		math.exp(-(standardDelay*w)**2/2) 
			for w in (2*math.pi*wavenumber*spectralRange/nBins for wavenumber in xrange(nBins) )
		])
	#computedSpectrum[abs(computedSpectrum)<.5e-2]=.5e-2
	computedInverseAudio = numpy.fft.irfft(1/computedSpectrum)

	bmaudio.saveWave("lala.wav",computedInverseAudio,44100)

	fixedSpectrum = numpy.array([
		math.sqrt(w)/10000+math.exp(-(standardDelay*w)**2/2)
			for w in (2*math.pi*wavenumber*spectralRange/nBins for wavenumber in xrange(nBins) ) ])

	import Gnuplot
	gp=Gnuplot.Gnuplot(persist=1)
	gp('set data style lines')
	gp("set log y ")
	gp.plot(
		numpy.zeros(nBins),
		abs(simulatedSpectrum),
		abs(computedSpectrum)[abs(computedSpectrum)>1e-4],
		abs(fixedSpectrum),
		abs(1/simulatedSpectrum),
		abs(1/computedSpectrum)[abs(computedSpectrum)>1e-4],
		abs(1/fixedSpectrum),
#		[bin.real for bin in spectrum],
#		[bin.imag for bin in spectrum],
	)
	gp.hardcopy(filename="IncoherenceSimulation.png",terminal="png") 


if not profile :
	main()
else: 
	import hotshot
	prof = hotshot.Profile("profile.valgrind")
	prof.runcall(main)
	prof.close()


import sympy

t=sympy.Symbol("t")
mu=sympy.Symbol("mu")
sigma=sympy.Symbol("s")

normal=1./(sigma*sympy.sqrt(2*sympy.pi))*sympy.exp(-(t/2/sigma)**2)
print "normal", normal
sympy.pretty_print(normal)
print sympy.printing.latex(normal)

sine=sympy.exp(sympy.I * 2*sympy.pi*mu*t/nBins)
print "sine", sine
sympy.pretty_print(sine)

print "product", normal*sine
sympy.pretty_print(normal*sine)

integral = sympy.integrate(normal*sine,t)
print "integral", integral

sympy.Plot(integral.subs(sigma,standardDelay*1))

