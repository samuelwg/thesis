#!/usr/bin/env python

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *

import headDistortion

def polarPatternX(azimuth, order) :
	return [ math.cos(order * azimuth) ]
def polarPatternY(azimuth, order) :
	return [ math.sin(order * azimuth) ]

nBins = nBins/16
spectrumBins = nBins/2+1
spectralRange=22050 # TODO: to debug
w=numpy.linspace(0, spectralRange*2*numpy.pi, spectrumBins, endpoint=True)

import sympy

N=720
colors="rgbcymb"

def adjustAngleForHeadDelay(azimuthRadians) :
	# moving the degrees to the -90,90 interval
	return azimuthRadians
	while azimuthRadians<  -math.pi : azimuthRadians += 2*math.pi
	while azimuthRadians>= +math.pi : azimuthRadians -= 2*math.pi
	# mirroring as it were in the first and fourth quadrant
	if azimuthRadians > math.pi/2 : azimuthRadians = math.pi-azimuthRadians
	if azimuthRadians < -math.pi/2 : azimuthRadians = -math.pi-azimuthRadians
	return azimuthRadians

def headDelaySimplified(azimuthRadians) :
	"""
	Returns the delay of a wave from a source to the ear considering the head difraction.
	Delay is computed relative to the delay to reach the head center without head.
	Formula has been simplified considering that head radius (R) << distance to source (r).
	Azimuth of the source is defined between -pi and pi. 0 at front, 90 at left (listening ear).
	"""
#	azimuthRadians = adjustAngleForHeadDelay(azimuthRadians)
	if azimuthRadians > 0 : return -R/c*math.sin(azimuthRadians)
	return -R/c*numpy.arcsin(numpy.sin(azimuthRadians))

def headlessDelaySimplified(azimuthRadians) :
	"""
	Returns the delay of a wave from a source to the ear position if there were no head.
	See headDelaySimplified for assumptions.
	"""
#	azimuthRadians = adjustAngleForHeadDelay(azimuthRadians)
	return -R/c*math.sin(azimuthRadians)

fullAzimuths = numpy.linspace(-math.pi, math.pi, num=N, endpoint=False) # All azimuths
directAzimuths = numpy.linspace(0, math.pi, num=N/2, endpoint=False) # Azimuths arround the listening ear
shadowAzimuths = numpy.linspace(-math.pi, 0, num=N/2, endpoint=False) # Azimuths arround the opposite ear

if True :
	print """Polar pattern of the delay of an spherical head:
	Checking the direct and shadow parts of the polar delay pattern.
	Should be near to a cardioid with minima at 90 (direct sound to the left ear).
	For the direct sound, in red it matches the cardioid.
	For the shadowed sound, it is a linear which follows the tangent of the
	cardioid and ends at (1+pi/2)*R/c.
	The polar graph is ofseted by R/c so it can be shown without negative parts.
	"""

	for label, line, azimuthSet in [
		('direct', 'r+', directAzimuths),
		('shadow', 'gx', shadowAzimuths),
		('full', '', fullAzimuths),
		] :
		pylab.polar(azimuthSet, R/c + numpy.array([headDelaySimplified(a) for a in azimuthSet]), line, label=label)
	pylab.polar(fullAzimuths, R/c + numpy.array([headlessDelaySimplified(a) for a in fullAzimuths]), line, label="headless")
	pylab.legend()
	pylab.show()


if True :
	print """Multiple angle polar patterns. Used by spherical harmonics.
	Absolute value displayed. Alternate lobules switch the sign.
	Components    0 degrees       90 degrees
	  X even       Symmetry        Symmetry
	  X odd:       Symmetry      Antisymmetry
	  Y even:    Antisymmetry    Antisymmetry
	  Y odd:     Antisymmetry      Symmetry
	Considering the integral for each quarter of the inner product of such 
	amplitude patterns with the delay pattern, which is symmetric around
	azimuth 90 gives (Qn integral in quarter n):
	- 90 Antisym. patterns (Xodd,Yeven) are cancelled (same delay opposite amplitudes)
	- 90 Symmetric patterns (Xeven,Yodd) Q1=Q2 and Q3=Q4
	If delay pattern is also antisymmetric at 0 (headless displaced case). This means:
	- X even: Q1=C(Q4) Q2=C(Q3) Q = Q1+Q2+Q3+Q4 = 4*Re(Q1)  (C: complementary)
	- Y odd:  Q1=-C(Q4) Q2=-C(Q3) Q = Q1+Q2+Q3+Q4= 4i*Im(Q1)
	"""
	def plotMultipleAnglePolarPattern(title, azimuths, polarPattern, ordersToShow) :
		pylab.title(title)
		for order in ordersToShow :
			pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
			pylab.polar(azimuths, numpy.abs(pattern), label=str(order))
		pylab.show()
	plotMultipleAnglePolarPattern("X even components (symmetry at Y)", fullAzimuths, polarPatternX, xrange(0,7,2))
	plotMultipleAnglePolarPattern("X odd components (anti-symmetry at Y)", fullAzimuths, polarPatternX, xrange(1,7,2)) # all zero
	plotMultipleAnglePolarPattern("Y even components (symmetry at Y)", fullAzimuths, polarPatternY, xrange(2,7,2)) # all zero
	plotMultipleAnglePolarPattern("Y odd components (anti-symmetry at Y)", fullAzimuths, polarPatternY, xrange(1,7,2))

if True :
	print """Headless displaced case:
	Without head diffraction just the displacement: Similar case which matches bessel function.
	X components:
	- Odd orders are zero,
	- Even orders is a 2*pi*BesselJ(order)
	Y components:
	- Even orders are zero,
	- Odd orders is i*2*pi*BesselJ(order)
	Because the symmetries and antisymmetries, this means:
	- X even:  Re(Q1) = Re(Q2) = BesselJ(order)*pi/2
	- Y odd:   Im(Q1) = Im(Q2) = BesselJ(order)*pi/2
	Q1 and Q2 matches with the direct sound part of our problem.
	So which are the other complex component of this part?
	"""
	def plotDisplacedResponseVsBesselFunctions(azimuths, delays, ordersToShow, polarPattern) :
		for order in ordersToShow :
			print order
			ocolor=colors[order%len(colors)]
			pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
			spectrum = headDistortion.delayAmplitudeToSpectrum(delays, pattern, w)
			besselN = numpy.array([float(sympy.mpmath.besselj(order, aw*R/c)) for aw in w ])
			pylab.plot(numpy.real(spectrum), ocolor+'-', label="R"+str(order))
			pylab.plot(numpy.imag(spectrum), ocolor+'--', label="I"+str(order))
			pylab.plot(besselN+0.01, ocolor+':', label="B"+str(order))
		pylab.ylim(-1,1)
		pylab.grid(.3)
		pylab.legend()
		pylab.show()
	delays = numpy.array([headlessDelaySimplified(a) for a in fullAzimuths])
	ordersToShow = xrange(0,5,1)
	plotDisplacedResponseVsBesselFunctions(fullAzimuths, delays, ordersToShow, polarPatternX)
	plotDisplacedResponseVsBesselFunctions(fullAzimuths, delays, ordersToShow, polarPatternY)

def fact2(x) :
	"Does not work for negatives, but enough"
	from operator import mul
	return reduce(mul, xrange(x,0,-2),1)


ordersToShow = xrange(0,16,1)

if True :
	"Checking struve symmetry ratio"
	orders = xrange(0,5)
	for order in orders :
		ocolor=colors[order%len(colors)]
		struveHdiff = lambda n,x : sympy.mpmath.struveh(n,x) - sympy.mpmath.struveh(-n,x)
		from operator import mul
		hypothesis = lambda n,x : \
			sum((
				2/numpy.pi * 
				( fact2(2*m+1)* x**(m-1) - (-1)**(n) * fact2(2*m-3) * x**(-m+1) ) *
				reduce(mul,((2*k-2)/x for k in xrange(n,m,-1)), 1)
				for m in xrange(n)
				))
		targetData = numpy.array([struveHdiff(order,x) for x in numpy.arange(0,100,1) ])
		target2Data = numpy.array([sympy.mpmath.struveh(order,x) for x in numpy.arange(0,100,1) ])
		target3Data = numpy.array([sympy.mpmath.struveh(-order,x) for x in numpy.arange(0,100,1) ])
		hypothesisData = numpy.array([hypothesis(order,x) for x in numpy.arange(0,100,1) ])
#		pylab.plot(hypothesisData, ocolor+"-", label="H"+str(order))
#		pylab.plot(targetData, ocolor+':', label="T"+str(order)) # should be zero
#		pylab.plot(target2Data, ocolor+'--', label="+"+str(order)) # should be zero
		pylab.plot(target3Data, ocolor+'', label="-"+str(order)) # should be zero
	pylab.legend()
	pylab.show()

if True :
	print """
	Checking the formulas for Rn.
	x=w*R/c
	1/pi int(sin(x sin(a) - n a), a,0,pi) =~
		1/pi sum(sin(x sin(a) - n a), a,0,pi,pi/N)/N*pi =
		sum(sin(x sin(a) - n a), a,0,pi,pi/N)/N 

	(-1)^n H_n(x) + 
	"""
	off=.1
	ordersToShow = xrange(0,15)
	x = w*(R/c) 
	for order in ordersToShow :
		print order
		ocolor=colors[order%len(colors)]
		target = sum([
			numpy.sin(x*numpy.sin(a)-order*a)
			for a in directAzimuths ])/len(directAzimuths)
		phase = (-1) ** order
		struveH = phase*numpy.array([float(sympy.mpmath.struveh(-order, ax)) for ax in x ])
		polyN = numpy.array([ 
			2/numpy.pi*
			sum((
				fact2(2*(order-m)-1)/fact2(2*m-1)*pow(ax,2*m-order-1)
				for m in xrange(1,order//2+1) ))
			for ax in x ])

		print order
		hypothesis = struveH + polyN

		remainder = hypothesis - target
#		pylab.plot(abs(polyN[1:]), ocolor+"-", label="Poly"+str(order))
#		pylab.plot(abs(struveH[1:]), ocolor+":*", label="Struve"+str(order))
		pylab.plot(remainder[1:], ocolor+"-", label="Remainder"+str(order)) # should be zero
		pylab.plot(target[1:], ocolor+'-x', label="Target"+str(order))
		pylab.plot(hypothesis[1:]+off, ocolor+'+', label="Hypotesis"+str(order))
#		pylab.plot(20*numpy.log10(numpy.abs(hypothesis)), ocolor, label=str(order))
	pylab.ylim(-4,4)
	pylab.grid(.3)
	pylab.legend()
	pylab.show()
	


ordersToShow = xrange(10,15,1)
if False :
	print """Direct sound component:
	To solve it we are missing Im(Q1)=Im(Q2) for Xeven and Re(Q1)=Re(Q2) for Yodd.
	Let x=wt/c
	For even orders:
		int(t,0,pi) cos(nt) exp(ixsin(t))
		= int(t,0,pi) cos(nt) exp(ixsin(t)) - i int(t,0,pi) sin(nt) exp(ixsin(t))
		= int(t,0,pi) (cos(nt)-isin(nt)) exp(ixsin(t))
		= int(t,0,pi) exp(-int) exp(ixsin(t))
	For odd orders:
		int(t,0,pi) sin(nt) exp(ixsin(t))
		= int(t,0,pi) sin(nt) exp(ixsin(t)) + i int(t,0,pi) cos(nt) exp(ixsin(t))
		= int(t,0,pi) (icos(nt)+sin(nt)) exp(ixsin(t))
		= i int(t,0,pi) (cos(nt)-isin(nt)) exp(ixsin(t))
		= i int(t,0,pi) exp(-int) exp(ixsin(t))
	Which is the same expression than for even orders but multiplied by i.
	So:
		int(t,0,pi) exp(-int) exp(ixsin(t))
		= int(t,-pi,pi) exp(-int) exp(ixsin(t)) - int(t,-pi,0) exp(-int) exp(ixsin(t))
		= 2pi Jn(x) - int(t,-pi,0) exp(-int) exp(ixsin(t))
		= 2pi Jn(x) + int(t,0,-pi) exp(-int) exp(ixsin(t))
		= 2pi Jn(x) - int(t,0,pi) exp(int) exp(-ixsin(t))
		= 2pi Jn(x) - int(t,0,pi) exp(int-ixsin(t))
		= 2pi Jn(x) - int(t,0,pi) cos(nt-xsin(t))+isin(nt-xsin(t))
		= 2pi Jn(x) - int(t,0,pi) cos(nt-xsin(t)) - i int(t,0,pi) sin(nt-xsin(t))
		= 2pi Jn(x) - pi Jn(x) - i int(t,0,pi) sin(nt-xsin(t))
		= pi Jn(x) - i int(t,0,pi) sin(nt-xsin(t))
		= pi Jn(x) + i int(t,0,pi) sin(xsin(t)-nt)

	By defining Rn(x) := 1/pi int(t,0,pi) sin(xsin(t)-nt)  (which is real)
		int(t,0,pi) exp(-int) exp(ixsin(t)) = pi Jn(x) + i pi Rn(x)
	so, for even orders:
		int(t,0,pi) cos(nt) exp(ixsin(t)) = pi (Jn(x) + i Rn(x))
	and for odd orders:
		int(t,0,pi) sin(nt) exp(ixsin(t)) = i pi (Jn(x) + i Rn(x))

	So for a given order n the only component is:
		i^((1-(-)^n)/2) pi (J_n(x)+iR_n(x)), where
		R_n(x) = 1/pi * \int_0^pi \sin(x \sin t - n t) dt

	The plot shows the difference between the numerical integration and 
	the analytic expression of the direct sound component. It should be zero.
	"""
	azimuths = directAzimuths
	delays = numpy.array([headDelaySimplified(a) for a in azimuths])
	off=0
#	ordersToShow = xrange(0,6,1)
	for order in ordersToShow :
		print order
		ocolor=colors[order%len(colors)]
		oddOrder = order&1

		polarPattern = polarPatternY if oddOrder else polarPatternX
		pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
		spectrum = headDistortion.delayAmplitudeToSpectrum(delays, pattern, w)

		phase = 1j**((1-(-1)**order)/2) # 1, j, 1, j...
		phase = 1j if order&1 else 1
		besselN = numpy.array([
			float(sympy.mpmath.besselj(order, aw*R/c))
			for aw in w ])
		otherN = sum([
			numpy.sin(w*R/c*numpy.sin(a)-order*a)
			for a in directAzimuths ])/len(directAzimuths) # per and divided by pi
		hypothesis = phase*(besselN + 1j*otherN)

		remainder = spectrum - hypothesis

		pylab.plot(hypothesis+.01, label="H"+str(order))
#		pylab.plot(numpy.real(remainder)+off, ocolor+':', label="R"+str(order)) # should be zero
#		pylab.plot(numpy.imag(remainder)+off, ocolor+'--', label="I"+str(order)) # should be zero
#		pylab.plot(besselN+off, label="B"+str(order)) ; off=off+.01
#		pylab.plot((-1)**order*besselN+off, label="B"+str(order)) ; off=off+.01
		pylab.plot(20*numpy.log10(numpy.abs(spectrum)), ocolor, label=str(order))
#	pylab.ylim(-1,1)
	pylab.grid(.3)
	pylab.legend()
	pylab.show()


if False :
	print """Shadowed sound component:
	For even orders:
		int(t,-pi,0) cos(nt) exp(ixasin(sin(t)))
		= int(t,-pi,0) cos(nt) exp(ixasin(sin(t))) - i int(t,0,pi) sin(nt) exp(ixasin(sin(t)))
		= int(t,-pi,0) (cos(nt)-isin(nt)) exp(ixasin(sin(t)))
		= int(t,-pi,0) exp(-int) exp(ixasin(sin(t)))
		= I_n(t)
	For odd orders:
		int(t,-pi,0) sin(nt) exp(ixasin(sin(t)))
		= int(t,-pi,0) sin(nt) exp(ixasin(sin(t))) + i int(t,0,pi) cos(nt) exp(ixasin(sin(t)))
		= int(t,-pi,0) (sin(nt)+icos(nt)) exp(ixasin(sin(t)))
		= i int(t,-pi,0) (-isin(nt)+cos(nt)) exp(ixasin(sin(t)))
		= i int(t,-pi,0) exp(-int) exp(ixasin(sin(t)))
		= i I_n(t)

	I_n(x) = int(t,-pi,0) exp(-int+ixasin(sin(t)))
	"""
	azimuths = shadowAzimuths
	delays = numpy.array([headDelaySimplified(a) for a in azimuths])
	off=0.01
#	ordersToShow = xrange(0,4,1)
	for order in ordersToShow :
		print order
		ocolor=colors[order%len(colors)]
		oddOrder = order&1

		polarPattern = polarPatternY if oddOrder else polarPatternX
		pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
		spectrum = headDistortion.delayAmplitudeToSpectrum(delays, pattern, w)

		azimuths0PiHalf = numpy.linspace(0,numpy.pi/2,N/4)
		azimuthsPiHalfPi = numpy.linspace(numpy.pi/2,numpy.pi,N/4)

		phase = 1j if order&1 else 1.
		# = 1/pi int(t,-pi,0) exp(-int+ixasin(sin(t)))
		# = 1/pi int(t,0,pi) exp(int-ixasin(sint))
		"""
		hypothesis = (phase/len(directAzimuths))*numpy.array([
			sum((
				numpy.exp(1j*order*a -1j*Rwc*numpy.arcsin(numpy.sin(a)))
				for a in directAzimuths
			)) for Rwc in (R/c*ws for ws in w) ]) 
		"""
		# = 1/pi int(t,0,pi/2) exp(int-ixt) + 1/pi int(t,pi/2,pi) exp(int-ix(pi-t))
		# 	T=pi-t
		# 	t=pi-T
		# 	dT=-dt
		# = 1/pi int(t,0,pi/2) exp(int-ixt) + 1/pi int(t,0,pi/2) exp(in(pi-t)-ixt)
		# = 1/pi int(t,0,pi/2) [ exp(int-ixt) + exp(in(pi-t)-ixt) ]
		"""
		hypothesis = (phase/len(azimuths0PiHalf)/2)*numpy.array([
			sum((
				numpy.exp(1j*order*a -1j*Rwc*a)
				+ numpy.exp(1j*order*(numpy.pi-a)-1j*Rwc*a)
				for a in azimuths0PiHalf
			)) for Rwc in (R/c*ws for ws in w) ])
		"""
		# = 1/pi int(t,0,pi/2) [ exp(int-ixt) + exp(inpi) exp(-int-ixt) ]
		# = 1/pi int(t,0,pi/2) [ exp(i(n-x)t) + exp(inpi) exp(-i(n+x)t) ]
		"""
		hypothesis = (phase/len(azimuths0PiHalf)/2)*numpy.array([
			sum((
				numpy.exp(1j*(order-Rwc)*a)
				+ numpy.exp(1j*order*numpy.pi) * numpy.exp(-1j*(order+Rwc)*a)
				for a in azimuths0PiHalf
			)) for Rwc in (R/c*ws for ws in w) ])
		"""
		# = 1/pi [0:pi/2] [ exp(i(n-x)t) / (i(n-x)) - exp(inpi) exp(-i(n+x)t) / (i(n+x)) ]
		# = 1/pi [0:pi/2] -i [ exp(i(n-x)t) / (n-x) - exp(inpi) exp(-i(n+x)t) / (n+x) ]
		"""
		hypothesis = (phase * -1j/numpy.pi )*(
			numpy.array([
				(
					+ numpy.exp(1j*(order-Rwc)*numpy.pi/2) / (order-Rwc)
					- numpy.exp(1j*order*numpy.pi) * numpy.exp(-1j*(order+Rwc)*numpy.pi/2) / (order+Rwc)
					- numpy.exp(1j*(order-Rwc)*0) / (order-Rwc)
					+ numpy.exp(1j*order*numpy.pi) * numpy.exp(-1j*(order+Rwc)*0) / (order+Rwc)
				) 
			for Rwc in R/c*w ])
			)
		"""
		# = -i/pi  [ ( exp(i(n-x)pi/2)-1 ) / (n-x) - exp(inpi) (exp(-i(n+x)pi/2)-1) / (n+x) ]
		# = -i/pi [ ( exp(i(n-x)pi/2)-1 ) / (n-x) - (exp(inpi-i(n+x)pi/2)-exp(inpi)) / (n+x) ]
		# = -i/pi [ ( exp(i(n-x)pi/2)-1 ) / (n-x) - (exp(inpi-inpi/2-ixpi/2)-exp(inpi)) / (n+x) ]
		# = -i/pi [ ( exp(i(n-x)pi/2)-1 ) / (n-x) - (exp(inpi/2-ixpi/2)-exp(inpi)) / (n+x) ]
		# = -i/pi [ ( exp(i(n-x)pi/2)-1 ) / (n-x) - (exp(i(n-x)pi/2)-exp(inpi)) / (n+x) ]
		"""
		hypothesis = (phase * -1j/numpy.pi )*(
			numpy.array([
				(
					+ numpy.exp(1j*(order-Rwc)*numpy.pi/2) / (order-Rwc)
					- 1. / (order-Rwc)
					- numpy.exp(1j*(order-Rwc)*numpy.pi/2) / (order+Rwc)
					+ numpy.exp(1j*order*numpy.pi) / (order+Rwc)
				) 
			for Rwc in R/c*w ])
			)
		"""
		# = -i/pi/(n*n-x*x) [ ( exp(i(n-x)pi/2)-1 ) * (n+x) - (exp(i(n-x)pi/2)-exp(inpi))* (n-x) ]
		# = -i/pi/(n*n-x*x) [ -(n+x) + (n-x)exp(inpi) + exp(i(n-x)pi/2) * 2x ]
		hypothesis = (phase * -1j/numpy.pi )*(
			numpy.array([
				(
					+ 2*Rwc* numpy.exp(1j*(order-Rwc)*numpy.pi/2)
					- (order+Rwc)
					+ (order-Rwc) * numpy.exp(1j*order*numpy.pi)
				) / (order*order-Rwc*Rwc)
			for Rwc in R/c*w ])
			)
		remainder = spectrum - hypothesis
#		pylab.plot(numpy.real(remainder), colors[order]+'--', label="R"+str(order))
#		pylab.plot(numpy.imag(remainder), colors[order]+':', label="I"+str(order))
#		pylab.plot(numpy.real(hypothesis)+off, colors[order]+'-', label="hR"+str(order))
#		pylab.plot(numpy.imag(hypothesis)+off, colors[order]+'.', label="hI"+str(order))
		pylab.plot(20*numpy.log10(numpy.abs(spectrum)), ocolor, label=str(order))
#	pylab.ylim(-30,1)
#	pylab.grid(.3)
#	pylab.legend()
#	pylab.show()

def analyticalComponent(order,w,t0) :
		phase = 1j if order&1 else 1.
		besselN = numpy.array([
			float(sympy.mpmath.besselj(order, wi*R/c))
			for wi in w ])
		otherN = sum([
			numpy.sin(w*R/c*numpy.sin(a)-order*a)
			for a in directAzimuths ])/len(directAzimuths)
		result = (phase * -1j/numpy.pi )* numpy.array([
				(
					+ 2*Rwc* numpy.exp(1j*(order-Rwc)*numpy.pi/2)
					- (order+Rwc)
					+ (order-Rwc) * numpy.exp(1j*order*numpy.pi)
				) / (order*order-Rwc*Rwc)
			for Rwc in t0*w ]) + phase*(besselN + 1j*otherN)
		if order is 0 : 
			result[0] = phase * (-1j + besselN[0] + 1j*otherN[0])
			print result[0], besselN[0], 1j*otherN[0]
#		result[numpy.isnan(result)] = numpy.pi #2.1
		return result

def numericalComponent(order, w, delays) :
		polarPattern = polarPatternY if order&1 else polarPatternX
		pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
		return headDistortion.delayAmplitudeToSpectrum(delays, pattern, w)*2
		# the 2 factor is added because we should multiply by the variable range (2pi)
		# while in the previous half cases we should have multiplied by pi.
		# We ignore the pi factor in all of them but keep the 2


if False :
	print """Whole thing in analytic
	"""
	azimuths = fullAzimuths
	delays = numpy.array([headDelaySimplified(a) for a in azimuths])
	off=0.01
#	ordersToShow = xrange(0,4,1)
	for order in ordersToShow :
		print order
		ocolor = colors[order%len(colors)]
		spectrum = numericalComponent(order, w, delays)
		hypothesis = analyticalComponent(order, w, R/c)

		remainder = spectrum - hypothesis
#		pylab.plot(numpy.real(remainder), ocolor+'--', label="RR"+str(order))
#		pylab.plot(numpy.imag(remainder), ocolor+':', label="RI"+str(order))
#		pylab.plot(numpy.real(spectrum), ocolor+'--', label="SR"+str(order))
#		pylab.plot(numpy.imag(spectrum), ocolor+':', label="SI"+str(order))
#		pylab.plot(numpy.real(hypothesis)+off, ocolor+'-', label="HR"+str(order))
#		pylab.plot(numpy.imag(hypothesis)+off, ocolor+'-.', label="HI"+str(order))
		pylab.plot(20*numpy.log10(numpy.abs(spectrum)), ocolor+'--', label=str(order))
#		pylab.plot(20*numpy.log10(numpy.abs(hypothesis)), ocolor+':', label=str(order))
	pylab.ylim(-1,1)
	pylab.ylim(-50,5)
	pylab.grid(.3)
	pylab.legend()
	pylab.show()

if False :
	"""Comparing the ifft of the analytic with the analitical time domain formula:
	+ 1/2pi int(a,-pi,pi) exp(-ina) delta(t-delay(a))
	:: split in two regions
	= 1/2pi *
		+ int(a,-pi,0) exp(-ina) delta(t-delay(a))
		+ int(a,0,pi) exp(-ina) delta(t-delay(a))
	:: different delay() for each region
	= 1/2pi *
		+ int(a,-pi,0) exp(-ina) delta(t+t0asin(sin(a)))
		+ int(a,0,pi)  exp(-ina) delta(t+t0sin(a))
	:: a'=-a and invert integral direction
	= 1/2pi *
		+ int(a,0,pi) exp(ina)  delta(t-t0asin(sin(a)))
		+ int(a,0,pi) exp(-ina) delta(t+t0sin(a))
	:: split in quarters
	= 1/2pi *
		+ int(a,0,pi/2)  exp(ina) delta(t-t0asin(sin(a)))
		+ int(a,pi/2,pi) exp(ina) delta(t-t0asin(sin(a)))
		+ int(a,0,pi/2)  exp(-ina) delta(t+t0sin(a))
		+ int(a,pi/2,pi) exp(-ina) delta(t+t0sin(a))
	:: a'=pi-a  a=pi-a' da = -da'
	= 1/2pi *
		+ int(a,0,pi/2) exp(ina) delta(t-t0asin(sin(a)))
		+ int(a,0,pi/2) exp(in(pi-a)) delta(t-t0asin(sin(pi-a)))
		+ int(a,0,pi/2) exp(-ina) delta(t+t0sin(a))
		+ int(a,0,pi/2) exp(-in(pi-a)) delta(t+t0sin(pi-a))
	:: sin pi-a = sin a
	= 1/2pi *
		+ int(a,0,pi/2) exp(ina) delta(t-t0asin(sin(a)))
		+ int(a,0,pi/2) exp(in(pi-a)) delta(t-t0asin(sin(a)))
		+ int(a,0,pi/2) exp(-ina) delta(t+t0sin(a))
		+ int(a,0,pi/2) exp(-in(pi-a)) delta(t+t0sin(a))
	:: asin(sin(a)) = a, en [0,pi/2]
	= 1/2pi *
		+ int(a,0,pi/2) exp(ina) delta(t-t0(a))
		+ int(a,0,pi/2) exp(in(pi-a)) delta(t-t0(a))
		+ int(a,0,pi/2) exp(-ina) delta(t+t0sin(a))
		+ int(a,0,pi/2) exp(-in(pi-a)) delta(t+t0sin(a))
	= 1/2pi *
		+ int(a,0,pi/2) exp(ina) delta(t-t0(a))
		+ int(a,0,pi/2) exp(-ina) delta(t-t0(a)) exp(inpi)
		+ int(a,0,pi/2) exp(-ina) delta(t+t0sin(a))
		+ int(a,0,pi/2) exp(ina) delta(t+t0sin(a)) exp(-inpi)
	:: v=sin(a) a=asin(v) da=dv/sqrt(1-a*a), v0=0, v1=1
	= 1/2pi *
		+ int(a,0,pi/2) exp(+ina) delta(t-t0(a))
		+ int(a,0,pi/2) exp(-ina) delta(t-t0(a)) exp(inpi)
		+ int(v,0,1) exp(-in asin(v)) delta(t+vt0) / sqrt(1-a*a)
		+ int(v,0,1) exp(+in asin(v)) delta(t+vt0) / sqrt(1-a*a) exp(-inpi)
	:: integrate (exp(1j*n*a)*deltafunctions.DiracDelta(t-a*t0), (a, 0, pi/2))
		= (Heaviside(t - pi*t0/2)- Heaviside(t)) * exp(I*n*t/t0)
	:: integrate (exp(-1j*n*a)*deltafunctions.DiracDelta(t-a*t0), (a, 0, pi/2))
		= (Heaviside(t - pi*t0/2)- Heaviside(t)) * exp(-I*n*t/t0)
	:: integrate (exp(1j*n*asin(a))*deltafunctions.DiracDelta(t-t0*a)/sqrt(1-a*a), (a, 0, 1))
		= (Heaviside(t - t0) - Heaviside(t)) * exp(I*n*asin(t/t0)) / (1 - t**2/t0**2)**(1/2)
	:: integrate (exp(-1j*n*asin(a))*deltafunctions.DiracDelta(t-t0*a)/sqrt(1-a*a), (a, 0, 1))
		= (Heaviside(t - t0) - Heaviside(t)) * exp(-I*n*asin(t/t0)) / (1 - t**2/t0**2)**(1/2)
	:: H(a,b) = Heaviside(t-a)-Heaviside(t-b) = Heaviside(t-a)*Heaviside(b-t)
	= 1/2pi *
		+ H(-t0 pi/2, 0) [
			+ exp(+int/t0)
			+ exp(-int/t0) exp(inpi)
			]
		+ H(0,t0) t0 / sqrt(t0t0 - tt) [
			+ exp(in asin(t/t0)) 
			+ exp(-in asin(t/t0)) exp(-inpi)
			]
	:: exp(iT)+exp(inpi)exp(-iT) = 2 cos(T) [n even) = 2i sin(T) [n odd]
	= 1/pi *
		+ H(-t0 pi/2, 0) [cos/isin](n t/t0)
		+ H(0,t0)   t0   [cos/isin](n asin(t/t0))  / sqrt(t0t0 - tt) 
	"""

	azimuths = fullAzimuths
	delays = numpy.array([headDelaySimplified(a) for a in azimuths])
	off=0.0
	ordersToShow = xrange(0,6,1)
	for order in ordersToShow :
		print order
		ocolor = colors[order%len(colors)]

		t0 = R/c
		analyticalFreqDomain = analyticalComponent(order, w, t0)
		reference = numpy.fft.irfft(analyticalFreqDomain)*18.
		
		phase = 1j if order&1 else 1.
		ir0 = [ t0/numpy.sqrt(t0*t0-t*t) * numpy.exp(1j*order*numpy.arcsin(t/t0))
			for t in numpy.linspace(-t0, 0, round(t0*samplingRate), endpoint=False) ]
		ir1 = [ numpy.exp(1j*order*t/t0)
			for t in numpy.linspace(0, t0*numpy.pi/2, round(t0*samplingRate/2*numpy.pi), endpoint=True) ]
		hypothesis = (phase*numpy.array(ir1+[0]*(nBins-len(ir0)-len(ir1))+ir0))

		reference = numpy.concatenate((reference[-round(t0*samplingRate)-1:-1],reference[0:round(t0*samplingRate/2*numpy.pi)+1] ))
		hypothesis = phase*numpy.array( ir0+ir1 )
		print len(hypothesis), len(reference), len(analyticalFreqDomain), len(w)
		remainder = reference #- hypothesis
		pylab.plot(numpy.real(remainder), ocolor+'--', label="RR"+str(order))
#		pylab.plot(numpy.imag(remainder), ocolor+':', label="RI"+str(order))
		pylab.plot(numpy.real(hypothesis)+off, ocolor+'-', label="HR"+str(order))
#		pylab.plot(numpy.imag(hypothesis)+off, ocolor+'-.', label="HI"+str(order))
#		pylab.plot(numpy.angle(spectrum), ocolor+'--', label=str(order))
#		pylab.plot(numpy.angle(hypothesis), ocolor+':', label=str(order))
#		pylab.plot(20*numpy.log10(numpy.abs(reference)), ocolor+'--', label=str(order))
#		pylab.plot(20*numpy.log10(numpy.abs(hypothesis)), ocolor+':', label=str(order))
#	pylab.ylim(-20,5)
	pylab.ylim(-2,2)
	pylab.grid(1)
	pylab.legend()
	pylab.show()



if False :
	"""Comparing the simplified delay with the one not considering R<<r.
	The remainder is multiplied by 100.
	It is a regular amplitude oscilation so that's why it is more significant
	on higher frequencies where the filter is lower, not because the absolute error increases.
	"""
	azimuths = fullAzimuths
	delays1 = numpy.array([headDelaySimplified(a) for a in azimuths])
	delays2 = numpy.array([bmaudio.sphericalHeadDelay(a,R,r)/c-r/c for a in azimuths])
	ordersToShow = xrange(0,1,1)
	for order in ordersToShow:
		print order
		ocolor=colors[order%len(colors)]
		oddOrder = order&1

		polarPattern = polarPatternY if oddOrder else polarPatternX
		pattern = [ polarPattern(azimuth, order) for azimuth in azimuths]
		spectrum1 = headDistortion.delayAmplitudeToSpectrum(delays1, pattern, w)
		spectrum2 = headDistortion.delayAmplitudeToSpectrum(delays2, pattern, w)

		remainder = spectrum1-spectrum2

		pylab.plot(100*numpy.real(remainder), ocolor+'--', label="RR"+str(order))
		pylab.plot(100*numpy.imag(remainder), ocolor+':', label="RI"+str(order))

		pylab.plot(20*numpy.log10(numpy.abs(spectrum1)), ocolor+"--", label="Aproximated "+str(order))
		pylab.plot(20*numpy.log10(numpy.abs(spectrum2)), ocolor+"-", label="Spherical "+str(order))
#		pylab.plot(numpy.angle(spectrum1), ocolor+"--", label="Aproximated "+str(order))
#		pylab.plot(numpy.angle(spectrum2), ocolor+"-", label="Spherical "+str(order))
	pylab.ylim(-numpy.pi,numpy.pi)
	pylab.ylim(-.1,.1)
	pylab.ylim(-40,5)
	pylab.grid(.3)
	pylab.legend()
	pylab.show()

if False :
	print """Plotting the filter along all the azimuths"""
	incidenceAzimuthsDegrees=numpy.linspace(0,180,40)
	incidenceAzimuths=numpy.linspace(0,numpy.pi,40)
	decoding=numpy.ones(3+1) # TODO: use common decodings, remove magic numbers
	spectrums = []
	delays = numpy.array([headlessDelaySimplified(a) for a in fullAzimuths])
	for incidenceIndex, incomingAzimuth in enumerate(incidenceAzimuths) :
		print incidenceIndex, incomingAzimuth
		spectrum = numpy.zeros(len(w))
		for order in xrange(3,4) :
			print "Order:",order
			oddOrder = order&1
			polarPattern = polarPatternY if oddOrder else polarPatternX
			pattern = [ polarPattern(incomingAzimuth, order)[0]* polarPattern(azimuth, order)[0] for azimuth in fullAzimuths]
			spectrum += headDistortion.delayAmplitudeToSpectrum(delays, pattern, w)*(decoding[order]/len(fullAzimuths)*2)
		spectrums.append(spectrum)
	print spectrums
	headDistortion.plotSpectrumsAsGreyMatrix(
		"Incidence dependant filter.eps", spectrums, spectralRange, incidenceAzimuthsDegrees )











