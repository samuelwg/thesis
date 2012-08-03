#!/usr/bin/env python
from __future__ import division
import glob, os, re, math, numpy, sys
import bmaudio
import pylab

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)

from parameters import *

def factorial(n) :
	return reduce(lambda x, y: x * y or 1, xrange(1,n+1), 1)

def semifactorial(n1, n2) :
	return float(reduce(lambda x, y: x * y or 1, xrange(n1,n2+1), 1))

def oldPolarPattern(a, NOrder) :
	return numpy.array([
		math.cos(n*a)
		for n in xrange(NOrder+1)
		])

def inphasePolarPattern(a, NOrder) :
	return numpy.array(
		[1] + [2*semifactorial(NOrder-n+1, NOrder) * math.cos(n*a) / semifactorial(NOrder+1,NOrder+n)
			for n in xrange(1,NOrder+1)] )

def basicPolarPattern(a, NOrder) :
	return numpy.array(
		[1] + [2*math.cos(n*a)
			for n in xrange(1,NOrder+1)] )

def maxrePolarPattern(a, NOrder) :
	return numpy.array(
		[1] + [2*math.cos(n*math.pi/(2*NOrder+2))*math.cos(n*a)
			for n in xrange(1,NOrder+1)] )

def nProduct(N,n) :
	import operator
	return reduce(operator.mul, xrange(N+1,N+n+1))

# TODO: unify decodingGains_XXX_nD and decodingsnDXXX

def decodingGains_MaxRv_2D(N) :
	return [1]*(N+1)
def decodingGains_MaxRv_3D(N) :
	return [1]*(N+1)
def decodingGains_InPhase_3D(N) :
	f = math.factorial
	return [ f(N) * f(N+1) / float( f(N+n+1) * f(N-n) ) for n in xrange(N+1) ]
	nP = nProduct
	return [ nP(N-n,n) / float( nP(N+1,n) ) for n in xrange(N+1) ]
def decodingGains_InPhase_2D(N) :
	f = math.factorial
	return [ f(N)**2 / float( f(N+n) * f(N-n) ) for n in xrange(N+1) ]
	nP = nProduct
	return [ nP(N-n,n) / float( nP(N,n) ) for n in xrange(N+1) ]

def decoding3dBasic(l) :
	return numpy.ones(l+1)

def decoding3dMaxRe(l) :
	if l is 0 : return numpy.array([1])
	import sympy
	x=sympy.Symbol("x")
	maxRoot = max((float(root) for root in sympy.solve(sympy.legendre(l+1,x),x)))
	return numpy.array([1.]+ [float(sympy.legendre(m, maxRoot))
		for m in xrange(1,l+1)])

def decoding3dInPhase(l) :
	gs = [1]+ [semifactorial(l-m+1,l)/semifactorial(l+2,l+m+1)
		for m in xrange(1,l+1)]
	return numpy.array(gs)

def decoding2dInPhase(l) :
	decoding = [1]+ [semifactorial(l-m+1,l)/semifactorial(l+1,l+m)
		for m in xrange(1,l+1)] # inphase 2D


w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

def spectrumDict(keys, spectrumBins) :
	return dict( 
		[ (key, numpy.zeros(spectrumBins,numpy.complex))
			for key in keys])

def sphericHeadMontecarlo(delayFunction, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow) :
	import random
	spectrums = spectrumDict(ordersToShow, spectrumBins)
	for i in xrange(speakers) :
		azimuthDegrees = random.random() * 360. # TODO: Check why 180 is not enough
		t = (delayFunction(math.radians(azimuthDegrees), R, r) - (r-R)) / c
		sinusoid = numpy.exp( 1j * t * w )
		a = math.radians(azimuthDegrees-referenceAzimuthDegrees)
		for order in ordersToShow :
			spectrums[order] += sum(polarPattern(a, order))/float(speakers) * sinusoid
	return spectrums

def delayBasedDecodingSimulation(delays, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow) :
	spectrums = spectrumDict(ordersToShow, spectrumBins)
	for azimuthDegrees in numpy.arange(speakers)*360/speakers :
		t = delays[azimuthDegrees]
		sinusoid = numpy.exp( 1j * t * w )
		a = math.radians(azimuthDegrees-referenceAzimuthDegrees)
		for order in ordersToShow :
			spectrums[order] += sum(polarPattern(a, order))/float(speakers) * sinusoid
	return spectrums

def sphericHeadSimulation(delayFunction, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow) :
	delays = dict([ 
		(azimuthDegrees, (delayFunction(math.radians(azimuthDegrees), R, r) - (r-R)) / c)
		for azimuthDegrees in numpy.arange(speakers)*360/speakers 
	])
	return delayBasedDecodingSimulation(delays, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow)

def hrtfDelaySimulation(hrtfDatabase, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow) :
	delays = bmaudio.azimuthDelaysFromHorizontalHrtf(hrtfDatabase)
	return delayBasedDecodingSimulation(delays, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow)

def hrtfFullSimulation(hrtfDatabase, referenceAzimuthDegrees, speakers, polarPattern, ordersToShow) :
	hrtfs = {}
	for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
		if elevationDegrees!= 0 : continue
		samplingRate, data = bmaudio.loadWave(response)
		hrtfs[azimuthDegrees] = numpy.fft.rfft(data,nBins)

	spectrums = spectrumDict(ordersToShow, spectrumBins)
	for azimuthDegrees in numpy.arange(speakers)*360/speakers :
		hrtf = hrtfs[azimuthDegrees]
		a = math.radians(azimuthDegrees-referenceAzimuthDegrees)
		for order in ordersToShow :
			spectrums[order] += sum(polarPattern(a,order))/float(speakers) * hrtf
	leftHrtf = hrtfs[referenceAzimuthDegrees]
	for order in ordersToShow :
		spectrums[order] /= leftHrtf
	return spectrums

def delayAmplitudeToSpectrum(delays, amplitudes, w) :
	"""
	Given an array of 'delays' (in radians) and their correspoinding 'amplitudes'
	constructs the spectrums of the resulting filter.
	w is an array containing the frequencies (rads/s) of each bin in the spectrum
	we want to get, and it is usually computed by:
	w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)
	"""
	return sum(a*numpy.exp(-1j*t*w) for t, a in zip(delays, amplitudes))/len(delays)

def plotSpectrumsAsGreyMatrix(name, spectrums, spectralRange, yaxis, yaxisLabel="Azimuth (degrees)", levels=30) :
	pylab.rcParams['figure.figsize']=(15,7)
	pylab.rcParams['figure.subplot.left'] = 0.05
	pylab.rcParams['figure.subplot.right'] = 0.90
	spectralBinPositions = numpy.arange(0,spectralRange, spectralRange/spectrumBins)
	pylab.contourf(spectralBinPositions, yaxis, spectrums,levels,cmap=pylab.cm.Greys_r, extend='both')
	pylab.colorbar().set_label("Magnitude (dB)")
	pylab.contour(spectralBinPositions, yaxis, spectrums,levels,cmap=pylab.cm.Greys)
	pylab.xlabel("Frequency (Hz)")
	pylab.ylabel(yaxisLabel)
	pylab.savefig(name, format="pdf")
	pylab.show()

if __name__ == "__main__" :
	sourceAzimuth = 90

	polarPattern = maxrePolarPattern
	plot = bmaudio.SpectrumDisplay()
	plot.inDb()
	#plot.showPhase()
#	plot.ylim(-20,5)
#	plot.flim(0,20000)
	ordersToShow = [1,3,5,7,9,15,17]
	ordersToShow = numpy.arange(7)+32
	N = 36
	for legend, filters in dict(
#		RealFull = hrtfFullSimulation(hrtfDatabase, sourceAzimuth, N, polarPattern, ordersToShow),
#		RealDelay = hrtfDelaySimulation(hrtfDatabase, sourceAzimuth, N, polarPattern, ordersToShow),
#		Spheric24 = sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, 24, polarPattern, ordersToShow),
#		Spheric36 = sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, 36, polarPattern, ordersToShow),
#		Spheric72 = sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, 72, polarPattern, ordersToShow),
		SphericAll = sphericHeadSimulation(bmaudio.sphericalHeadDelay, sourceAzimuth, 72, polarPattern, ordersToShow),
#		DelayAll = hrtfDelaySimulation(hrtfDatabase, sourceAzimuth, 72, polarPattern, ordersToShow),
	).items() :
		for order in ordersToShow :
			plot.addSpectrumData(
				filters[order], 22050, "%s o%i"%(legend,order))
	#plot.addSpectrumData(hrtfs[90], samplingRate/2, "90")
	plot.hardcopy(os.path.splitext(__file__)[0]+".pdf", "pdf")

def sphericalHead3dSH(order, samplingRate, nBins, R=R, c=c) :
	"""
	Returns the time domain 3D nth order component spherical harmonic
	of the HRTF of an spherical head at the left ear."""
	import sympy
	Pl = sympy.special.polynomials.legendre
	t0 = R/c
	duration = float(nBins)/samplingRate
	timeAxis = [x for x in numpy.linspace(-t0, duration-t0, nBins)]
#	timeAxis = [-t0 + float(i)/samplingRate for i in xrange(nBins)]
	return	numpy.sqrt(2*order+1)/2./t0 * numpy.array( [
			0. if t>math.pi/2*t0 else (
			Pl(order, (t/t0)) if (t<0) else
			numpy.cos(t/t0) * Pl(order, math.sin(t/t0)) 
			)
		for t in timeAxis ])

def sphericalHead3dSHSpectrum(order, spectralRange, nSpectrumBins, R=R, c=c) :
	from sympy import symbols, exp, sin, cos, I, simplify
	import operator
	import sympy
	t0=R/c

	w,t=symbols("w t")

	E=exp(-I*t*w)
	def N(n) : return sympy.Integer(n)
	def R(n) : return reduce(operator.mul, xrange(n,1,-2),sympy.Integer(1))
	def sum(v) : return reduce(operator.add, v, sympy.Integer(0))
	def f(n) : return reduce(operator.mul, xrange(n,0,-1), sympy.Integer(1))
	w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12 = [ w**(i+1) for i in xrange(12) ]
	t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12 = [ t**(i+1) for i in xrange(12) ]
	r1, r3, r5, r7, r9, r11, r13, r15, r17, r19, r21 = [R(i) for i in xrange(1,22,2) ]
	r0, r2, r4, r6, r8, r10, r12, r14, r16, r18, r20 = [R(i) for i in xrange(0,21,2) ]
	b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18 = [ sympy.Integer(2)**(i) for i in xrange(18+1) ]
	f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12 = [ f(i) for i in xrange(12+1) ]
	c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = [ cos(sympy.Integer(i)*t) for i in xrange(12+1) ]
	s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12 = [ sin(sympy.Integer(i)*t) for i in xrange(12+1) ]
	z0, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12 = [ (w**2-sympy.Integer(i)**2) for i in xrange(12+1) ]
	shadowCoeffs = [[
	 - N(1),
	],[
	 + r1/b1,
	],[
	 + r1/r2/b2,
	 + r3/f2/b2,
	],[
	 - r3/f3/b2,
	 - r5/f3/b3,
	],[
	 + r3/N(3)/b6,
	 - r5/b7,
	 - r7/f4/b4,
	],[
	 - r5/N(3)/b8,
	 + r7/N(3)/N(5)/b6,
	 + r9/f5/b5,
	],[
	 + r5/N(3)/b10,
	 - r7/N(5)/b10,
	 + r9/N(3)/N(3)/b10,
	 + r11/f6/b6,
	],[
	 - r7/f5/b7,
	 + r9/f6/b6,
	 - r11*N(3)/f7/b6,
	 - r13/f7/b7,
	],[
	 + r7/f3/b13,
	 - r9/f5/b10,
	 + r11/f4/N(7)/N(3)/b10,
	 - r13/f6/b6/b5,
	 - r15/f8/b8,
	],[
	 - r9/N(9)/b15,
	 + r11/N(9)/N(7)/N(5)/b12,
	 - r13/f7/b12*N(3),
	 + r15/f9/b6,
	 + r17/f9/b9,
	],[
	 + r9/N(9)/N(5)*N(7)/b17,
	 - r11/N(7)/N(3)/b17,
	 + r13/N(9)/N(7)/b18,
	 - r15/f6/b6/N(9)/b8,
	 + r17/f8/b8/N(5)/b3,
	 + r19/f10/b10,
	]]

	formulaShadow = sum( [
		coefficient * (
			i*sin(i*t)-I*w*cos(i*t)
			if i&1  else
			i*cos(i*t)+I*w*sin(i*t)
			)
			/(w*w-i*i)
		for i, coefficient in zip(xrange(1+(order&1),order+2,2), shadowCoeffs[order])
		])

	formulaDirect = I**(1-order) * sum([
		I**wi
		* w**wi 
		* sum([
			(-1)**( wi//2-N(i) )
			* simplify(
				R(2*order -1-wi+ti)
				/ R(wi-ti)
				/ factorial(ti)
			)
			* t**ti
			for i, ti in zip(xrange(wi//2+1), xrange(wi&1, wi+1, 2))
			])
		for wi in xrange(order+1)
		])
	fullFormula = simplify(
		(2*order+1)**(1/N(2))/N(2) *( # TODO and still misses fft factor 1/2pi
			+ simplify(formulaShadow.subs(t,sympy.pi/N(2))*E.subs(t,sympy.pi/N(2)) - formulaShadow.subs(t,N(0)))
			+ simplify(formulaDirect.subs(t,N(0))/w**N(order+1) - formulaDirect.subs(t,N(-1))*E.subs(t,N(-1))/w**N(order+1))
		))

	ws=numpy.linspace(0,spectralRange*2*numpy.pi,nSpectrumBins,endpoint=False)
	return numpy.array([ complex(fullFormula.subs(w,wi*t0)) for wi in ws ])



