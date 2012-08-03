#!/usr/bin/env python

import glob, os, re, math, numpy
from bmaudio import *
import pylab
import random


print "Simulating horizontal sources an spheric head..."

N = 400000
r=16.0
R=0.075
a0 = math.acos(R/r)
nBins = 4*1024
samplingRate = 88200. # samples/s
c = 344.
d0 = math.sqrt(r*r-R*R)
tmax = (math.pi/2+1)*R/c

def analyticToni_Omni(i) :
	t = i/samplingRate
	if not i : t=.0000029
	if t>=tmax: return 0.
	result = 1./math.sqrt(2*R*t/c - t*t) if t<R/c else c/R
	return result/ math.pi/samplingRate

def analyticDavid_Omni(i) :
	t = i/samplingRate
	if not i : t=.0000029
	if t>=tmax: return 0.
	d = t*c-R+r
	result = 2.*(t*c-R+r)*c/math.sqrt(-(t*c+2*r)*(t*c-2*R+2*r)*(t*c)*(t*c-2*R)) if d<=math.sqrt(r*r-R*R) else c/R
	return result/math.pi/samplingRate

def analyticToni_Sin(i) :
	t = i/samplingRate
	if t>tmax: return 0.
	result = c/R *math.cos(1-t*c/R) if t>R/c else c/R
	return result/math.pi/samplingRate

def analyticDavid_Sin(i) :
	t = i/samplingRate
	if t>tmax: return 0.
	d = t*c-R+r
	if d<d0 :
		result = d/R/r
	else :
		result = math.sin( (d-d0)/R + a0) / R
	return c*result/math.pi/samplingRate

def analyticDavid_Cos(i) :
	t = i/samplingRate
	if not i : t=.0000029
	if t>tmax: return 0.
	d = t*c-R+r
	if d<d0 :
		rootContent = -(d-R-r)*(d-R+r)*(d+R-r)*(d+R+r)
		if (rootContent<=0) : return .0
		result = d * (R*R + r*r - d*d) / (R*r*math.sqrt(rootContent))
	else :
		result = math.cos((d-d0)/R+a0) /R
	return c*result/math.pi/samplingRate



simulated_Omni = numpy.zeros( nBins )
simulated_Cos = numpy.zeros( nBins )
simulated_Sin = numpy.zeros( nBins )
simulated_CosSin = numpy.zeros( nBins )
simulated_Sin2 = numpy.zeros( nBins )
for i in xrange(N) :
	i = random.random() * N # Montecarlo sampling (instead of equidistant)
	a = i * numpy.pi / N
	bin = (sphericalHeadDelay(a+math.pi/2, R, r) - (r-R)) * samplingRate / c
	simulated_Omni[bin] += 1./N
	simulated_Sin[bin] += math.sin(a)/N
	simulated_Cos[bin] += math.cos(a)/N
	simulated_CosSin[bin] += math.cos(a)*math.sin(a)/N
	simulated_Sin2[bin] += math.sin(a)**2/N

filterLeft = simulated_Omni + simulated_Cos
filterRight = simulated_Omni - simulated_Cos

filters = dict(
	[(f.func_name, numpy.array([f(i) for i in xrange(nBins)])) for f in [
		analyticDavid_Omni,
		analyticToni_Omni,
#		analyticDavid_Sin,
#		analyticToni_Sin,
#		analyticDavid_Cos,
]])

filters.update(
	simulated_Omni = simulated_Omni,
	simulated_Sin = simulated_Sin,
#	simulated_Cos = simulated_Cos,
#	simulated_CosSin = simulated_CosSin,
#	simulated_Sin2 = simulated_Sin2,
#	filterLeft = filterLeft,
#	filterRight = filterRight,
)


def dbFFT(data) :
	return 20*numpy.log10(numpy.abs(numpy.fft.rfft(data)))
def magFFT(data) :
	return numpy.abs(numpy.fft.rfft(data))


for label, data in filters.iteritems() :
	pylab.plot(data[:60], label=label)
pylab.legend()
pylab.show()


for label, data in filters.iteritems() :
	pylab.plot(dbFFT(data), label=label)
pylab.legend()
pylab.show()







