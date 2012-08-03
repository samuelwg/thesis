#!/usr/bin/python

import sympy
import numpy

azimuthSteps = 100
tSteps = 100

def function(t, azimuth) :
	from numpy import pi,sin,cos
	return (
		sin(4*2*pi*t) +
		sin(8*2*pi*t) *
		sin(4*2*pi*azimuth) +
		sin(1*2*pi*azimuth+4) +
		sin(2*2*pi*azimuth+2) +
		sin(3*2*pi*azimuth+3) +
		sin(9*2*pi*azimuth) +
		sin(11*2*pi*azimuth+pi) +
		sin(10*2*pi*azimuth) +
		sin(4*2*pi*azimuth) +
		sin(30*2*pi*azimuth+3) +
		sin(20*2*pi*azimuth+3) +
		sin(25*2*pi*azimuth+3) +
		sin(50*2*pi*t+3) +
		0
	)


result = numpy.zeros((tSteps, azimuthSteps))
for nt in xrange(tSteps) :
	for nazimuth in xrange(azimuthSteps) :
		result[nt, nazimuth] = function(float(nt)/tSteps, float(nazimuth)/azimuthSteps)



spectral = sum(numpy.fft.rfft(result[:][i]) for i in xrange(azimuthSteps))

import pylab

pylab.contour(result, 0,
	cmap=pylab.cm.RdBu)
pylab.contourf(result, 9,
	cmap=pylab.cm.RdBu,
#	levels = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5],
#	colors = ('w','k',)
	)
pylab.show()

pylab.plot(abs(spectral))
pylab.show()


