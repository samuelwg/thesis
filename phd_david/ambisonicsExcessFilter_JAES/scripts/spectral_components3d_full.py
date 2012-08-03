#!/usr/bin/env python
from __future__ import division
from spectral_components3d_common import *
from spectral_components3d_Pli import Pli
from spectral_components3d_direct import DirectDefiniteIntegral
from spectral_components3d_shadow import ShadowDefiniteIntegral

pi = sympy.pi

"""
The problem findng the spectral domain expression of Hl,0,+:
sqrt(2*l+1)/2/to * ( # TODO: maybe some 1/2pi is needed for fourier but this makes H0(0) be 0dB
	integrate(Pl(l,t/to) * exp(-I*w*t), (t,-to,0)) +
	integrate(cos(t/to) * Pl(l,sin(t/to)) * exp(-I*w*t), (t,0,to*pi/2)) 
	)
	Substitution:
		t = to * z
		dt = to * dz
		z0=-1 z1=0
		z0=0, z1= pi/2
sqrt(2*l+1)/2 * (
	+ integrate(Pl(l,z)*exp(-1j*w*to*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(l,sin(z))*exp(-1j*w*to*z), (z,0,pi/2)) # shadowed part
	)
	Substitution:
		y = -1j*w*to
sqrt(2*l+1)/2 * (
	+ integrate(Pl(l,z)*exp(y*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(l,sin(z))*exp(y*z), (z,0,pi/2) # shadowed part
	)
"""


maz.check("FullHl1", "takingformer, %i", [
	N(2*l+1)**(1/N(2))/N(2) * (ShadowDefiniteIntegral(l,w) + DirectDefiniteIntegral(l,w))
for l in xrange(10) ])

maz.check("FullHl1", "inlining former, %i", [
	N(2*l+1)**(1/N(2))/N(2) *
	(
		N(-1)**N(l)
		* N(2)
		*
		sum([
		sum([
			0 if r&1 or (l+1+n)&1 else
			I**N(r)
			* N(n) 
			* R(2*l-1-r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
			/ R(r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-sympy.pi/N(2)*I*w*to)
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			+ w*to  *(I**N(n)-I**N(-n)) / N(2)
		) / ((-I*w*to)**2+N(n*n))
		for n in xrange(1,l+1+1) ])
		+
		sum([
		(
			(0 if m&1 else 1)
			* I**N(m)
			* R(2*l-m-1)
			/ R(m)
			- exp(I*w*to) *
			sum([
				(0 if i&1 else 1)
				* I**N(i)
				* R(2*l-i-1)
				/ R(i)
				/ f(m-i)
			for i in xrange(m+1)])
		)
		/ I**N(l+m-1)
		/ (w*to)**N(l-m+1)
		for m in xrange(l+1) ])
	)
for l in xrange(8) ])

maz.check("FullHl1", "I^m phasors simplified, %i", [
	N(2*l+1)**(1/N(2))/N(2) *
	N(-1)**N(l) *
	(
		N(2)
		*
		sum([
		sum([
			0 if r&1 or (l+1+n)&1 else
			N(n) 
			* I**N(r)
			* R(2*l-1-r)
			/ R(r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-I*w*to*sympy.pi/N(2))
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			+ w*to  *(I**N(n)-I**N(-n)) / N(2)
		) / (+N(n)**2-(w*to)**2)
		for n in xrange(1,l+1+1) ])
		+
		sum([
		(
			(0 if m&1 else 1)
			* I**N(m)
			* R(2*l-1-m)
			/ R(m)
			- exp(I*w*to) *
			sum([
				(0 if i&1 else 1)
				* I**N(i)
				* R(2*l-1-i)
				/ R(i)
				/ f(m-i)
			for i in xrange(m+1)])
		)
		/ (-I*w*to)**N(l+1-m)
		for m in xrange(l+1) ])
	)
for l in xrange(8) ])

def Bcoeff(l,n) :
	return (
		(I**N(n)+I**N(-n))
		/ N(2)
		* R(2*l -1 -n)
		/ R(n)
		)

def Acoeff(l,n,r) :
	return (
			(I**N(l+1+n)+I**N(-l-1-n))**2
			/ N(4)
			* Bcoeff(l,r)
			/ R(l+1+n-r)
			/ R(l+1-n-r)
		)


maz.check("FullHl1", "Using B coefficients, %i", [
	N(2*l+1)**(1/N(2))/N(2) *
	N(-1)**N(l) *
	(
		N(2)
		*
		sum([
		sum([
			N(n) 
			* Acoeff(l,n,r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-I*w*to*sympy.pi/N(2))
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			+ w*to  *(I**N(n)-I**N(-n)) / N(2)
		) / (+N(n)**2-(w*to)**2)
		for n in xrange(1,l+1+1) ])
		+
		sum([
		(
			Bcoeff(l,m)
			- exp(I*w*to) *
			sum([
				Bcoeff(l,i)
				/ f(m-i)
			for i in xrange(m+1)])
		)
		/ (-I*w*to)**N(l+1-m)
		for m in xrange(l+1) ])
	)
for l in xrange(8) ])

maz.check("FullHl1", "Using B coefficients, %i", [
	N(2*l+1)**(1/N(2))/N(2) *
	N(-1)**N(l) *
	(
		N(2)
		*
		sum([
		sum([
			N(n) 
			* Acoeff(l,n,r)
		for r in xrange(0, l+1-n+1) ])
		* (
			- n     *(-1)**N(n) * exp(-I*w*to*sympy.pi/N(2))
			+ n     *(I**N(n)+I**N(-n)) / N(2)
			+ w*to  *(I**N(n)-I**N(-n)) / N(2)
		) / (+N(n)**2-(w*to)**2)
		for n in xrange(1,l+1+1) ])
		+
		sum([
		(
			Bcoeff(l,m)
			- exp(I*w*to) *
			sum([
				Bcoeff(l,i)
				/ f(m-i)
			for i in xrange(m+1)])
		)
		/ (-I*w*to)**N(l+1-m)
		for m in xrange(l+1) ])
	)
for l in xrange(8) ])


def debug(args) : print args; return N(1)

maz.check("Order 2", "integral as is, %i", [
	N(2*l+1)**(1/N(2))/N(2)/to * (
		(integrate(Pl(l,t/to) * sympy.exp(-I*w*t), (t,-to,0))) +
		(integrate(sympy.cos(t/to) * Pl(l,sympy.sin(t/to)) * sympy.exp(-I*w*t), (t,0,to*sympy.pi/2)))
	)
for l in xrange(0,2) ]+[
	N(2*l+1)**(1/N(2))/N(2) * (ShadowDefiniteIntegral(l,w) + DirectDefiniteIntegral(l,w))
for l in xrange(2,3) ])


print "Integrated"
for formula in maz.recall("Order 2") : print formula; print

maz.check("Order 2", "integral as is, %i", [
(
	+ I/w
	+ I*w*to**2/(
		+ 1
		- 2*to**2*w**2
		+ to**4*w**4
		)
	- I*to**4*w**3/(
		+ 1
		- 2*to**2*w**2
		+ to**4*w**4
		)
	- I*exp(I*to*w)/w
	- to*exp(pi*I*to*w/2)/(
		- exp(pi*I*to*w)
		- to**4*w**4*exp(pi*I*to*w)
		+ 2*to**2*w**2*exp(pi*I*to*w)
		)
	+ to**3*w**2*exp(pi*I*to*w/2)/(
		- exp(pi*I*to*w)
		- to**4*w**4*exp(pi*I*to*w)
		+ 2*to**2*w**2*exp(pi*I*to*w)
		)
)/(2*to)
,
N(3)**(1/N(2))*(
	-(exp(I*to*w)/w**2 - I*to*exp(I*to*w)/w)/to
	+ N(16)*I*to/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	+ I*to**5*w**4/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	- N(8)*I*to**3*w**2/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	+ 1/(to*w**2)
	+ N(16)*I*to*exp(pi*I*to*w)/(
		+ N(64)*I*exp(3*pi*I*to*w/2)
		- N(48)*I*to**2*w**2*exp(3*pi*I*to*w/2)
		+ N(12)*I*to**4*w**4*exp(3*pi*I*to*w/2)
		- N(1)*I*to**6*w**6*exp(3*pi*I*to*w/2)
		)
	+ I*to**5*w**4*exp(pi*I*to*w)/(
		+ N(64)*I*exp(3*pi*I*to*w/2)
		- N(48)*I*to**2*w**2*exp(3*pi*I*to*w/2)
		+ N(12)*I*to**4*w**4*exp(3*pi*I*to*w/2)
		- N(1)*I*to**6*w**6*exp(3*pi*I*to*w/2)
		)
	- N(8)*I*to**3*w**2*exp(pi*I*to*w)/(
		+ N(64)*I*exp(3*pi*I*to*w/2)
		- N(48)*I*to**2*w**2*exp(3*pi*I*to*w/2)
		+ N(12)*I*to**4*w**4*exp(3*pi*I*to*w/2)
		- N(1)*I*to**6*w**6*exp(3*pi*I*to*w/2)
		)
)/(2*to)
,
N(5)**(1/N(2))*(
	(
		- N(54)*I*to**0*w**0
		+ N(1)*I*to**4*w**4
		+ N(51)*I*to**2*w**2
		- N(2)*to**5*w**5*exp(-pi*I*to*w/2)
		- N(2)*I*to**6*w**6*exp(I*to*w)
		+ N(6)*to**5*w**5*exp(I*to*w)
		+ N(26)*I*to**4*w**4*exp(I*to*w)
		- N(60)*to**3*w**3*exp(I*to*w)
		- N(78)*I*to**2*w**2*exp(I*to*w)
		+ N(54)*to**1*w**1*exp(I*to*w)
		+ N(54)*I*to**0*w**0*exp(I*to*w)
	)/(
		+ N(2)*to**7*w**7
		- N(20)*to**5*w**5
		+ N(18)*to**3*w**3
	)
)/N(2)
])

maz.check("Order 2", "exp out of denominator, %i", [
(
	+ I/w
	- I*exp(I*to*w)/w
	- I*w*to**2/(
		- 1
		+ 2*to**2*w**2
		- to**4*w**4
		)
	+ I*to**4*w**3/(
		- 1
		+ 2*to**2*w**2
		- to**4*w**4
		)
	- to*exp(-pi*I*to*w/2)/(
		- 1
		+ 2*to**2*w**2
		- to**4*w**4
		)
	+ to**3*w**2*exp(-pi*I*to*w/2)/(
		- 1
		+ 2*to**2*w**2
		- to**4*w**4
		)
)/(2*to)
,
N(3)**(1/N(2))*(
	+ 1/(to*w**2)
	-(exp(I*to*w)/w**2 - I*to*exp(I*to*w)/w)/to
	+ N(16)*I*to/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	+ I*to**5*w**4/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	- N(8)*I*to**3*w**2/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	+ N(16)*I*to*exp(-pi*I*to*w/2)/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	+ I*to**5*w**4*exp(-pi*I*to*w/2)/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
	- N(8)*I*to**3*w**2*exp(-pi*I*to*w/2)/(
		+ N(64)*I
		- N(48)*I*to**2*w**2
		+ N(12)*I*to**4*w**4
		- N(1)*I*to**6*w**6
		)
)/(2*to)
,
N(5)**(1/N(2))*(

	(
		- N(54)*I*to**0*w**0
		+ N(1)*I*to**4*w**4
		+ N(51)*I*to**2*w**2
		- N(2)*to**5*w**5*exp(-pi*I*to*w/2)
		- N(2)*I*to**6*w**6*exp(I*to*w)
		+ N(6)*to**5*w**5*exp(I*to*w)
		+ N(26)*I*to**4*w**4*exp(I*to*w)
		- N(60)*to**3*w**3*exp(I*to*w)
		- N(78)*I*to**2*w**2*exp(I*to*w)
		+ N(54)*to**1*w**1*exp(I*to*w)
		+ N(54)*I*to**0*w**0*exp(I*to*w)
	)/(
		N(2)
		* (N(1)**2 - (w*to)**2)
		* (N(3)**2 - (w*to)**2)
		* (w*to)**3
	)
)/N(2)
])

maz.check("Order 2", "joining terms, %i", [
(
	+ I*(1-exp(I*to*w))/w**1/to**1
	+
		(I*w*to + exp(-pi*I*to*w/2))
		* (
			- 1
			+ to**2*w**2
		)/(
			- 1
			+ 2*to**2*w**2
			- to**4*w**4
			)
)/2
,
N(3)**(1/N(2))*(
	+ (1+(I*w*to - 1)*exp(I*to*w))/to**2/w**2
	+ 
		(1 + exp(-pi*I*to*w/2))
		*
		(
			+ N(16)
			- N(8)*to**2*w**2
			+ N(1)*to**4*w**4
		) /(
			+ N(64)
			- N(48)*to**2*w**2
			+ N(12)*to**4*w**4
			- N(1)*to**6*w**6
			)
)/2
,
N(5)**(1/N(2))*(

	+ (-I*w*to)**2 * exp(-pi*I*to*w/2)
		/(
		  (N(1)**2 - (w*to)**2)
		* (N(3)**2 - (w*to)**2)
		)
	- (
		+ N(1)*(-I*w*to)**2
		+ N(3)*(-I*w*to)**1
		+ N(3)*(-I*w*to)**0
	) * exp(I*to*w)
	/ (
		(-I*w*to)**3
	)
	+ (
		- N(1) * (-I*w*to)**2
		+ N(6)
	) / (
		(-I*w*to)**3
	) / N(2)

	+ (
		+ N(1) *(-I*w*to)
	)/(
		N(8)
		* (N(1)**2 - (w*to)**2)
	)
	+ (
		+ N(3) *(-I*w*to)
	)/(
		N(8)
		* (N(3)**2 - (w*to)**2)
	)
)/N(2)
])

maz.check("Order 2", "Numerator is denominator divided by (n**2-(w**n*to**n)) , %i", [
(
	+
		(1-exp(I*to*w))
		/ (-I*w*to)**1
	+
		(
			- N(1) * (-I*w*to)**1
			+ N(1) * (-I*w*to)**0 * exp(-pi*I*to*w/2)
		)
		/ (N(1) + (-I*w*to)**2)
)/2
,
N(3)**(1/N(2))*(
	+
		(-1+(1-I*w*to)*exp(I*to*w))
		/ (-I*w*to)**2
	+ 
		(
			+ N(1) * (-I*w*to)**0
			+ N(1) * (-I*w*to)**0 * exp(-pi*I*to*w/2)
		)
		/ (N(4) + (-I*w*to)**2)
)/2
,
N(5)**(1/N(2))*(

	+ (
		- N(1) * (-I*w*to)**2
		+ N(6)
	) / (
		(-I*w*to)**3
	) / N(2)
	+ (
		- N(1)*(-I*w*to)**2
		- N(3)*(-I*w*to)**1
		- N(3)*(-I*w*to)**0
	) * exp(I*to*w)
	/ (
		(-I*w*to)**3
	)

	+ (
		+ N(1) * (-I*w*to)**1
		+ N(1) * (-I*w*to)**2 * exp(-pi*I*to*w/2)
	)/(
		N(8)
		* (N(1)**2 + (-I*w*to)**2)
	)
	+ (
		+ N(3) * (-I*w*to)**1
		- N(1) * (-I*w*to)**2 * exp(-pi*I*to*w/2)
	)/(
		N(8)
		* (N(3)**2 + (-I*w*to)**2)
	)
)/N(2)
])

maz.check("Order 2", "Comparing concrete integral with general formula , %i", [
	N(2*l+1)**(1/N(2))/N(2) * (ShadowDefiniteIntegral(l,w) + DirectDefiniteIntegral(l,w))
for l in xrange(3) ])

print "Literal"
for formula in maz.recall("Order 2") : print formula; print

sys.exit(0)




for formula in maz.recall("FullHl1") : print formula; print

maz.check("FullHl1", "I^m phasors simplified, %i", [
	N(2*l+1)**(1/N(2))/N(2)/to * sympy.simplify(
		integrate(Pl(l,t/to) * sympy.exp(-I*w*t), (t,-to,0)) +
		integrate(sympy.cos(t/to) * Pl(l,sympy.sin(t/to)) * sympy.exp(-I*w*t), (t,0,to*sympy.pi/2))
	)
for l in xrange(2) ])

for formula in maz.recall("FullHl1") : print formula; print


sys.exit()


"""
# TODO: This is not something workable yet, although it is paperable
# By simplication of mathematica/sympy output


#target_shadow=[ integrate(Pl(sin(x))*cos(x)*exp(y*x),x)/exp(y*x) for l in xrange(10) ] # sympy cannot cope with it
Table[Integrate[Cos[t] LegendreP[l, Sin[t]] Exp[-I w t], t] /Exp(-I t w) , {l, 0, 10}]
{
(I w cos(t) - sin(t))/(-1 + w**2),
( 2 cos(2 t) + I w sin(2 t))/(2 (-4 + w**2)),
(1/( 8 (9 - 10 w**2 + w**4)))(-I w (-9 + w**2) cos(t) - 3 I w (-1 + w**2) cos(3 t) + 2 (-9 + 5 w**2 + 9 (-1 + w**2) cos(2 t)) sin( t)),
-(2 (-16 + w**2) cos(2 t) + 10 (-4 + w**2) cos( 4 t) + I w (-16 + w**2 + 5 (-4 + w**2) cos(2 t)) sin( 2 t))/(8 (64 - 20 w**2 + w**4)),
1/128 (-((2 I w cos(t))/(-1 + w**2)) + ( 15 I w cos(3 t))/(-9 + w**2) + ( 35 I w cos(5 t))/(-25 + w**2) + ( 2 sin(t))/(-1 + w**2) - (45 sin(3 t))/(-9 + w**2) - ( 175 sin(5 t))/(-25 + w**2)),
1/256 (-((10 cos(2 t))/(-4 + w**2)) + (112 cos(4 t))/(-16 + w**2) + ( 378 cos(6 t))/(-36 + w**2) - ( 5 I w sin(2 t))/(-4 + w**2) + ( 28 I w sin(4 t))/(-16 + w**2) + ( 63 I w sin(6 t))/(-36 + w**2)),
(1/1024)(-(( 5 I w cos(t))/(-1 + w**2)) + ( 21 I w cos(3 t))/(-9 + w**2) - ( 105 I w cos(5 t))/(-25 + w**2) - ( 231 I w cos(7 t))/(-49 + w**2) + ( 5 sin(t))/(-1 + w**2) - (63 sin(3 t))/(-9 + w**2) + ( 525 sin(5 t))/(-25 + w**2) + ( 1617 sin(7 t))/(-49 + w**2)),
(1/2048)(-(( 28 cos(2 t))/(-4 + w**2)) + (168 cos(4 t))/(-16 + w**2) - ( 1188 cos(6 t))/(-36 + w**2) - (3432 cos(8 t))/(-64 + w**2) - ( 14 I w sin(2 t))/(-4 + w**2) + ( 42 I w sin(4 t))/(-16 + w**2) - ( 198 I w sin(6 t))/(-36 + w**2) - ( 429 I w sin(8 t))/(-64 + w**2)),
(1/32768)(-(( 70 I w cos(t))/(-1 + w**2)) + ( 252 I w cos(3 t))/(-9 + w**2) - ( 660 I w cos(5 t))/(-25 + w**2) + ( 3003 I w cos(7 t))/(-49 + w**2) + ( 6435 I w cos(9 t))/(-81 + w**2) + ( 70 sin(t))/(-1 + w**2) - (756 sin(3 t))/(-9 + w**2) + ( 3300 sin(5 t))/(-25 + w**2) - (21021 sin(7 t))/(-49 + w**2) - ( 57915 sin(9 t))/(-81 + w**2)),
(1/65536)(-(( 420 cos(2 t))/(-4 + w**2)) + (2112 cos(4 t))/(-16 + w**2) - ( 7722 cos(6 t))/(-36 + w**2) + (45760 cos(8 t))/(-64 + w**2) + ( 121550 cos(10 t))/(-100 + w**2) - ( 210 I w sin(2 t))/(-4 + w**2) + ( 528 I w sin(4 t))/(-16 + w**2) - ( 1287 I w sin(6 t))/(-36 + w**2) + ( 5720 I w sin(8 t))/(-64 + w**2) + ( 12155 I w sin(10 t))/(-100 + w**2)),
(1/262144)(-(( 294 I w cos(t))/(-1 + w**2)) + ( 990 I w cos(3 t))/(-9 + w**2) - ( 2145 I w cos(5 t))/(-25 + w**2) + ( 5005 I w cos(7 t))/(-49 + w**2) - ( 21879 I w cos(9 t))/(-81 + w**2) - ( 46189 I w cos(11 t))/(-121 + w**2) + ( 294 sin(t))/(-1 + w**2) - (2970 sin(3 t))/(-9 + w**2) + ( 10725 sin(5 t))/(-25 + w**2) - (35035 sin(7 t))/(-49 + w**2) + ( 196911 sin(9 t))/(-81 + w**2) + (508079 sin(11 t))/(-121 + w**2))
}
"""

target_shadow = [ # minimum translation to python (syntax, cos and sin)
(I*w*c1 - s1)/(-1 + w**2),
( 2*c2 + I*w*s2)/(2*(-4 + w**2)),
(1/8/(9 - 10*w**2 + w**4))*(-I*w*(-9 + w**2)*c1 - 3*I*w*(-1 + w**2)*c3 + (-9+ w**2)*s1 + 9 * s3 *(-1 + w**2)),
-(2 *(-16 + w**2) *cos(2* t) + 10* (-4 + w**2)* cos( 4* t) + I* w* (-16 + w**2 + 5 *(-4 + w**2) *cos(2* t)) *sin( 2* t))/(8 *(64 - 20* w**2 + w**4)),
1/128*(-((2*I*w*c1)/(-1 + w**2)) + ( 15*I*w*c3)/(-9 + w**2) + ( 35*I*w*c5)/(-25 + w**2) + ( 2*s1)/(-1 + w**2) - (45*s3)/(-9 + w**2) - ( 175*s5)/(-25 + w**2)),
1/256*(-((10*c2)/(-4 + w**2)) + (112*c4)/(-16 + w**2) + ( 378*c6)/(-36 + w**2) - ( 5*I*w*s2)/(-4 + w**2) + ( 28*I*w*s4)/(-16 + w**2) + ( 63*I*w*s6)/(-36 + w**2)),
(1/1024)*(-(( 5*I*w*c1)/(-1 + w**2)) + ( 21*I*w*c3)/(-9 + w**2) - ( 105*I*w*c5)/(-25 + w**2) - ( 231*I*w*c7)/(-49 + w**2) + ( 5*s1)/(-1 + w**2) - (63*s3)/(-9 + w**2) + ( 525*s5)/(-25 + w**2) + ( 1617*s7)/(-49 + w**2)),
(1/2048)*(-(( 28*c2)/(-4 + w**2)) + (168*c4)/(-16 + w**2) - ( 1188*c6)/(-36 + w**2) - (3432*c8)/(-64 + w**2) - ( 14*I*w*s2)/(-4 + w**2) + ( 42*I*w*s4)/(-16 + w**2) - ( 198*I*w*s6)/(-36 + w**2) - ( 429*I*w*s8)/(-64 + w**2)),
(1/32768)*(-(( 70*I*w*c1)/(-1 + w**2)) + ( 252*I*w*c3)/(-9 + w**2) - ( 660*I*w*c5)/(-25 + w**2) + ( 3003*I*w*c7)/(-49 + w**2) + ( 6435*I*w*c9)/(-81 + w**2) + ( 70*s1)/(-1 + w**2) - (756*s3)/(-9 + w**2) + ( 3300*s5)/(-25 + w**2) - (21021*s7)/(-49 + w**2) - ( 57915*s9)/(-81 + w**2)),
(1/65536)*(-(( 420*c2)/(-4 + w**2)) + (2112*c4)/(-16 + w**2) - ( 7722*c6)/(-36 + w**2) + (45760*c8)/(-64 + w**2) + ( 121550*c10)/(-100 + w**2) - ( 210*I*w*s2)/(-4 + w**2) + ( 528*I*w*s4)/(-16 + w**2) - ( 1287*I*w*s6)/(-36 + w**2) + ( 5720*I*w*s8)/(-64 + w**2) + ( 12155*I*w*s10)/(-100 + w**2)),
(1/262144)*(-(( 294*I*w*c1)/(-1 + w**2)) + ( 990*I*w*c3)/(-9 + w**2) - ( 2145*I*w*c5)/(-25 + w**2) + ( 5005*I*w*c7)/(-49 + w**2) - ( 21879*I*w*c9)/(-81 + w**2) - ( 46189*I*w*c11)/(-121 + w**2) + ( 294*s1)/(-1 + w**2) - (2970*s3)/(-9 + w**2) + ( 10725*s5)/(-25 + w**2) - (35035*s7)/(-49 + w**2) + ( 196911*s9)/(-81 + w**2) + (508079*s11)/(-121 + w**2))
]

simplifiedShadow = [  # after doing many simplifications
 - 1* (1*s1 - I*w*c1)/z1
,
 + r1/b1* (2*c2 + I*w*s2)/z2
,
 + r1/r2/b2* (1*s1 - I*w*c1)/z1
 + r3/f2/b2* (3*s3 - I*w*c3)/z3
,
 - r3/f3/b2* (2*c2 + I*w*s2)/z2
 - r5/f3/b3* (4*c4 + I*w*2*s2*c2)/z4 # 2*s2*c2 == s4, but sympy does not see it
,
 + r3/3/b6*   (1*s1 - I*w*c1)/z1
 - r5/b7*     (3*s3 - I*w*c3)/z3
 - r7/f4/b4*  (5*s5 - I*w*c5)/z5
,
 - r5/3/b8*   (2*c2 + I*w*s2)/z2
 + r7/3/5/b6* (4*c4 + I*w*s4)/z4
 + r9/f5/b5*  (6*c6 + I*w*s6)/z6
,
 + r5/3/b10*   (1*s1 - I*w*c1)/z1
 - r7/5/b10*   (3*s3 - I*w*c3)/z3
 + r9/3/3/b10* (5*s5 - I*w*c5)/z5
 + r11/f6/b6*  (7*s7 - I*w*c7)/z7
,
 - r7/f5/b7*    (2*c2 + I*w*s2)/z2
 + r9/f6/b6*    (4*c4 + I*w*s4)/z4
 - r11*3/f7/b6* (6*c6 + I*w*s6)/z6
 - r13/f7/b7*   (8*c8 + I*w*s8)/z8
,
 + r7/f3/b13*     (1*s1 - I*w*c1)/z1
 - r9/f5/b10*     (3*s3 - I*w*c3)/z3
 + r11/f4/7/3/b10*(5*s5 - I*w*c5)/z5
 - r13/f6/b6/b5*  (7*s7 - I*w*c7)/z7
 - r15/f8/b8*     (9*s9/z9 - I*w*c9/z9) # TODO: common factor z9
,
 - r9/9/b15*      (2*c2 + I*w*s2)/z2
 + r11/9/7/5/b12* (4*c4 + I*w*s4)/z4
 - r13/f7/b12*3*  (6*c6 + I*w*s6)/z6
 + r15/f9/b6*     (8*c8 + I*w*s8)/z8
 + r17/f9/b9*     (10*c10 + I*w*s10)/z10
,
 + r9/9/5*7/b17*   (1*s1 - I*w*c1)/z1
 - r11/7/3/b17*    (3*s3 - I*w*c3)/z3
 + r13/9/7/b18*    (5*s5 - I*w*c5)/z5
 - r15/f6/b6/9/b8* (7*s7 - I*w*c7)/z7
 + r17/f8/b8/5/b3* (9*s9/z9 - I*w*c9/z9) # TODO: common factor z9
 + r19/f10/b10*    (11*s11/z11 - I*w*c11/z11) # TODO common factor z11
]
simplifiedShadow2 = [ # clearing the things sympy does not see as equal
 - 1* (1*s1 - I*w*c1)/z1
,
 + r1/b1* (2*c2 + I*w*s2)/z2
,
 + r1/r2/b2* (1*s1 - I*w*c1)/z1
 + r3/f2/b2* (3*s3 - I*w*c3)/z3
,
 - r3/f3/b2* (2*c2 + I*w*s2)/z2
 - r5/f3/b3* (4*c4 + I*w*s4)/z4
,
 + r3/3/b6*   (1*s1 - I*w*c1)/z1
 - r5/b7*     (3*s3 - I*w*c3)/z3
 - r7/f4/b4*  (5*s5 - I*w*c5)/z5
,
 - r5/3/b8*   (2*c2 + I*w*s2)/z2
 + r7/3/5/b6* (4*c4 + I*w*s4)/z4
 + r9/f5/b5*  (6*c6 + I*w*s6)/z6
,
 + r5/3/b10*   (1*s1 - I*w*c1)/z1
 - r7/5/b10*   (3*s3 - I*w*c3)/z3
 + r9/3/3/b10* (5*s5 - I*w*c5)/z5
 + r11/f6/b6*  (7*s7 - I*w*c7)/z7
,
 - r7/f5/b7*    (2*c2 + I*w*s2)/z2
 + r9/f6/b6*    (4*c4 + I*w*s4)/z4
 - r11*3/f7/b6* (6*c6 + I*w*s6)/z6
 - r13/f7/b7*   (8*c8 + I*w*s8)/z8
,
 + r7/f3/b13*     (1*s1 - I*w*c1)/z1
 - r9/f5/b10*     (3*s3 - I*w*c3)/z3
 + r11/f4/7/3/b10*(5*s5 - I*w*c5)/z5
 - r13/f6/b6/b5*  (7*s7 - I*w*c7)/z7
 - r15/f8/b8*     (9*s9 - I*w*c9)/z9 # TODO: common factor z9
,
 - r9/9/b15*        (2*c2 + I*w*s2)/z2
 + r11/9/7/5/b12*      (4*c4 + I*w*s4)/z4
 - r13/f7/b12*3*   (6*c6 + I*w*s6)/z6
 + r15/f9/b6*   (8*c8 + I*w*s8)/z8
 + r17/f9/b9*     (10*c10 + I*w*s10)/z10
,
 + r9/9/5*7/b17*   (1*s1 - I*w*c1)/z1
 - r11/7/3/b17*    (3*s3 - I*w*c3)/z3
 + r13/9/7/b18*    (5*s5 - I*w*c5)/z5
 - r15/f6/b6/9/b8* (7*s7 - I*w*c7)/z7
 + r17/f8/b8/5/b3* (9*s9 - I*w*c9)/z9 # TODO: common factor z9
 + r19/f10/b10*    (11*s11 - I*w*c11)/z11 # TODO common factor z11
]
shadowCoeffs = [[
 - 1,
],[
 + r1/b1,
],[
 + r1/r2/b2,
 + r3/f2/b2,
],[
 - r3/f3/b2,
 - r5/f3/b3,
],[
 + r3/3/b6,
 - r5/b7,
 - r7/f4/b4,
],[
 - r5/3/b8,
 + r7/3/5/b6,
 + r9/f5/b5,
],[
 + r5/3/b10,
 - r7/5/b10,
 + r9/3/3/b10,
 + r11/f6/b6,
],[
 - r7/f5/b7,
 + r9/f6/b6,
 - r11*3/f7/b6,
 - r13/f7/b7,
],[
 + r7/f3/b13,
 - r9/f5/b10,
 + r11/f4/7/3/b10,
 - r13/f6/b6/b5,
 - r15/f8/b8,
],[
 - 9*7*5/3/2**15,
 + 11*9/3/2**12,
 - 13*11*9/2**16,
 + 15*13*11/3/2**13,
 + 17*15*13*11/3/2**16,
],[
 + r9/9/5*7/b17,
 - r11/7/3/b17,
 + r13/9/7/b18,
 - r15/f6/b6/9/b8,
 + r17/f8/b8/5/b3,
 + r19/f10/b10,
]]


formulaShadow = [
	sum( [
		coefficient * (
			i*sin(i*t)-I*w*cos(i*t)
			if i&1  else
			i*cos(i*t)+I*w*sin(i*t)
			)
			/(w*w-i*i)
		for i, coefficient in zip(xrange(1+(n&1),n+2,2), coefficientsN)
		])
	for n,coefficientsN in enumerate(shadowCoeffs)
]


def shadowCoefficientsHipothesis(n,i) :
	return 0

for n, orderCoeficients in enumerate(shadowCoeffs) :
	print n, [simplify(sympy.sympify(x-shadowCoefficientsHipothesis(n,i))) for i,x in enumerate(reversed(orderCoeficients))]

# Test: should be [0,0,0...]
print "diff shadow:", [simplify(target_shadow[i] - simplifiedShadow[i]) for i in xrange(11)]

# Test: should be [0,0,0...]
print "diff shadow formula:", [simplify(simplifiedShadow2[i] - formulaShadow[i]) for i in xrange(11)]


"""
################################################################
# Direct part integral
################################################################
"""

#target_direct = [ integrate(E*special.polynomials.legendre(i,t),t).doit()*w**(i+1) / E for i in xrange(12) ]
#print target_direct
# precomputed
target_direct = [I, w**2*(exp(-I*t*w)/w**2 + I*t*exp(-I*t*w)/w)*exp(I*t*w), w**3*(-3*I*exp(-I*t*w)/w**3 + 3*t*exp(-I*t*w)/w**2 - I*exp(-I*t*w)/(2*w) + 3*I*t**2*exp(-I*t*w)/(2*w))*exp(I*t*w), w**4*(-15*exp(-I*t*w)/w**4 - 3*exp(-I*t*w)/(2*w**2) + 15*t**2*exp(-I*t*w)/(2*w**2) - 15*I*t*exp(-I*t*w)/w**3 - 3*I*t*exp(-I*t*w)/(2*w) + 5*I*t**3*exp(-I*t*w)/(2*w))*exp(I*t*w), w**5*(-105*t*exp(-I*t*w)/w**4 + 105*I*exp(-I*t*w)/w**5 - 15*t*exp(-I*t*w)/(2*w**2) + 3*I*exp(-I*t*w)/(8*w) + 15*I*exp(-I*t*w)/(2*w**3) + 35*t**3*exp(-I*t*w)/(2*w**2) - 105*I*t**2*exp(-I*t*w)/(2*w**3) - 15*I*t**2*exp(-I*t*w)/(4*w) + 35*I*t**4*exp(-I*t*w)/(8*w))*exp(I*t*w), w**6*(945*exp(-I*t*w)/w**6 + 15*exp(-I*t*w)/(8*w**2) + 105*exp(-I*t*w)/(2*w**4) - 945*t**2*exp(-I*t*w)/(2*w**4) - 105*t**2*exp(-I*t*w)/(4*w**2) + 315*t**4*exp(-I*t*w)/(8*w**2) + 945*I*t*exp(-I*t*w)/w**5 - 315*I*t**3*exp(-I*t*w)/(2*w**3) - 35*I*t**3*exp(-I*t*w)/(4*w) + 15*I*t*exp(-I*t*w)/(8*w) + 63*I*t**5*exp(-I*t*w)/(8*w) + 105*I*t*exp(-I*t*w)/(2*w**3))*exp(I*t*w), w**7*(-10395*I*exp(-I*t*w)/w**7 + 10395*t*exp(-I*t*w)/w**6 - 3465*t**3*exp(-I*t*w)/(2*w**4) - 945*I*exp(-I*t*w)/(2*w**5) - 315*t**3*exp(-I*t*w)/(4*w**2) - 105*I*exp(-I*t*w)/(8*w**3) - 5*I*exp(-I*t*w)/(16*w) + 105*t*exp(-I*t*w)/(8*w**2) + 693*t**5*exp(-I*t*w)/(8*w**2) + 945*t*exp(-I*t*w)/(2*w**4) - 3465*I*t**4*exp(-I*t*w)/(8*w**3) - 315*I*t**4*exp(-I*t*w)/(16*w) + 105*I*t**2*exp(-I*t*w)/(16*w) + 231*I*t**6*exp(-I*t*w)/(16*w) + 945*I*t**2*exp(-I*t*w)/(4*w**3) + 10395*I*t**2*exp(-I*t*w)/(2*w**5))*exp(I*t*w), w**8*(-135135*exp(-I*t*w)/w**8 - 10395*exp(-I*t*w)/(2*w**6) - 945*exp(-I*t*w)/(8*w**4) - 35*exp(-I*t*w)/(16*w**2) - 45045*t**4*exp(-I*t*w)/(8*w**4) - 3465*t**4*exp(-I*t*w)/(16*w**2) + 945*t**2*exp(-I*t*w)/(16*w**2) + 3003*t**6*exp(-I*t*w)/(16*w**2) + 10395*t**2*exp(-I*t*w)/(4*w**4) + 135135*t**2*exp(-I*t*w)/(2*w**6) - 135135*I*t*exp(-I*t*w)/w**7 - 10395*I*t*exp(-I*t*w)/(2*w**5) - 9009*I*t**5*exp(-I*t*w)/(8*w**3) - 945*I*t*exp(-I*t*w)/(8*w**3) - 693*I*t**5*exp(-I*t*w)/(16*w) - 35*I*t*exp(-I*t*w)/(16*w) + 315*I*t**3*exp(-I*t*w)/(16*w) + 429*I*t**7*exp(-I*t*w)/(16*w) + 3465*I*t**3*exp(-I*t*w)/(4*w**3) + 45045*I*t**3*exp(-I*t*w)/(2*w**5))*exp(I*t*w), w**9*(-2027025*t*exp(-I*t*w)/w**8 + 2027025*I*exp(-I*t*w)/w**9 - 135135*t*exp(-I*t*w)/(2*w**6) - 135135*t**5*exp(-I*t*w)/(8*w**4) - 10395*t*exp(-I*t*w)/(8*w**4) - 9009*t**5*exp(-I*t*w)/(16*w**2) - 315*t*exp(-I*t*w)/(16*w**2) + 35*I*exp(-I*t*w)/(128*w) + 315*I*exp(-I*t*w)/(16*w**3) + 3465*t**3*exp(-I*t*w)/(16*w**2) + 6435*t**7*exp(-I*t*w)/(16*w**2) + 10395*I*exp(-I*t*w)/(8*w**5) + 45045*t**3*exp(-I*t*w)/(4*w**4) + 135135*I*exp(-I*t*w)/(2*w**7) + 675675*t**3*exp(-I*t*w)/(2*w**6) - 2027025*I*t**2*exp(-I*t*w)/(2*w**7) - 135135*I*t**2*exp(-I*t*w)/(4*w**5) - 45045*I*t**6*exp(-I*t*w)/(16*w**3) - 10395*I*t**2*exp(-I*t*w)/(16*w**3) - 3003*I*t**6*exp(-I*t*w)/(32*w) - 315*I*t**2*exp(-I*t*w)/(32*w) + 3465*I*t**4*exp(-I*t*w)/(64*w) + 6435*I*t**8*exp(-I*t*w)/(128*w) + 45045*I*t**4*exp(-I*t*w)/(16*w**3) + 675675*I*t**4*exp(-I*t*w)/(8*w**5))*exp(I*t*w), w**10*(34459425*exp(-I*t*w)/w**10 + 315*exp(-I*t*w)/(128*w**2) + 3465*exp(-I*t*w)/(16*w**4) + 135135*exp(-I*t*w)/(8*w**6) + 2027025*exp(-I*t*w)/(2*w**8) - 34459425*t**2*exp(-I*t*w)/(2*w**8) - 2027025*t**2*exp(-I*t*w)/(4*w**6) - 765765*t**6*exp(-I*t*w)/(16*w**4) - 135135*t**2*exp(-I*t*w)/(16*w**4) - 45045*t**6*exp(-I*t*w)/(32*w**2) - 3465*t**2*exp(-I*t*w)/(32*w**2) + 45045*t**4*exp(-I*t*w)/(64*w**2) + 109395*t**8*exp(-I*t*w)/(128*w**2) + 675675*t**4*exp(-I*t*w)/(16*w**4) + 11486475*t**4*exp(-I*t*w)/(8*w**6) + 34459425*I*t*exp(-I*t*w)/w**9 - 11486475*I*t**3*exp(-I*t*w)/(2*w**7) - 675675*I*t**3*exp(-I*t*w)/(4*w**5) - 109395*I*t**7*exp(-I*t*w)/(16*w**3) - 45045*I*t**3*exp(-I*t*w)/(16*w**3) - 6435*I*t**7*exp(-I*t*w)/(32*w) - 1155*I*t**3*exp(-I*t*w)/(32*w) + 315*I*t*exp(-I*t*w)/(128*w) + 3465*I*t*exp(-I*t*w)/(16*w**3) + 9009*I*t**5*exp(-I*t*w)/(64*w) + 12155*I*t**9*exp(-I*t*w)/(128*w) + 135135*I*t*exp(-I*t*w)/(8*w**5) + 135135*I*t**5*exp(-I*t*w)/(16*w**3) + 2027025*I*t*exp(-I*t*w)/(2*w**7) + 2297295*I*t**5*exp(-I*t*w)/(8*w**5))*exp(I*t*w), w**11*(-654729075*I*exp(-I*t*w)/w**11 + 654729075*t*exp(-I*t*w)/w**10 - 218243025*t**3*exp(-I*t*w)/(2*w**8) - 34459425*I*exp(-I*t*w)/(2*w**9) - 11486475*t**3*exp(-I*t*w)/(4*w**6) - 2078505*t**7*exp(-I*t*w)/(16*w**4) - 2027025*I*exp(-I*t*w)/(8*w**7) - 675675*t**3*exp(-I*t*w)/(16*w**4) - 109395*t**7*exp(-I*t*w)/(32*w**2) - 45045*I*exp(-I*t*w)/(16*w**5) - 15015*t**3*exp(-I*t*w)/(32*w**2) - 3465*I*exp(-I*t*w)/(128*w**3) - 63*I*exp(-I*t*w)/(256*w) + 3465*t*exp(-I*t*w)/(128*w**2) + 45045*t*exp(-I*t*w)/(16*w**4) + 135135*t**5*exp(-I*t*w)/(64*w**2) + 230945*t**9*exp(-I*t*w)/(128*w**2) + 2027025*t*exp(-I*t*w)/(8*w**6) + 2297295*t**5*exp(-I*t*w)/(16*w**4) + 34459425*t*exp(-I*t*w)/(2*w**8) + 43648605*t**5*exp(-I*t*w)/(8*w**6) - 218243025*I*t**4*exp(-I*t*w)/(8*w**7) - 11486475*I*t**4*exp(-I*t*w)/(16*w**5) - 2078505*I*t**8*exp(-I*t*w)/(128*w**3) - 675675*I*t**4*exp(-I*t*w)/(64*w**3) - 109395*I*t**8*exp(-I*t*w)/(256*w) - 15015*I*t**4*exp(-I*t*w)/(128*w) + 3465*I*t**2*exp(-I*t*w)/(256*w) + 45045*I*t**2*exp(-I*t*w)/(32*w**3) + 45045*I*t**6*exp(-I*t*w)/(128*w) + 46189*I*t**10*exp(-I*t*w)/(256*w) + 765765*I*t**6*exp(-I*t*w)/(32*w**3) + 2027025*I*t**2*exp(-I*t*w)/(16*w**5) + 14549535*I*t**6*exp(-I*t*w)/(16*w**5) + 34459425*I*t**2*exp(-I*t*w)/(4*w**7) + 654729075*I*t**2*exp(-I*t*w)/(2*w**9))*exp(I*t*w), w**12*(-13749310575*exp(-I*t*w)/w**12 - 654729075*exp(-I*t*w)/(2*w**10) - 34459425*exp(-I*t*w)/(8*w**8) - 675675*exp(-I*t*w)/(16*w**6) - 45045*exp(-I*t*w)/(128*w**4) - 693*exp(-I*t*w)/(256*w**2) - 4583103525*t**4*exp(-I*t*w)/(8*w**8) - 218243025*t**4*exp(-I*t*w)/(16*w**6) - 43648605*t**8*exp(-I*t*w)/(128*w**4) - 11486475*t**4*exp(-I*t*w)/(64*w**4) - 2078505*t**8*exp(-I*t*w)/(256*w**2) - 225225*t**4*exp(-I*t*w)/(128*w**2) + 45045*t**2*exp(-I*t*w)/(256*w**2) + 675675*t**2*exp(-I*t*w)/(32*w**4) + 765765*t**6*exp(-I*t*w)/(128*w**2) + 969969*t**10*exp(-I*t*w)/(256*w**2) + 14549535*t**6*exp(-I*t*w)/(32*w**4) + 34459425*t**2*exp(-I*t*w)/(16*w**6) + 305540235*t**6*exp(-I*t*w)/(16*w**6) + 654729075*t**2*exp(-I*t*w)/(4*w**8) + 13749310575*t**2*exp(-I*t*w)/(2*w**10) - 13749310575*I*t*exp(-I*t*w)/w**11 - 916620705*I*t**5*exp(-I*t*w)/(8*w**7) - 654729075*I*t*exp(-I*t*w)/(2*w**9) - 43648605*I*t**5*exp(-I*t*w)/(16*w**5) - 34459425*I*t*exp(-I*t*w)/(8*w**7) - 4849845*I*t**9*exp(-I*t*w)/(128*w**3) - 2297295*I*t**5*exp(-I*t*w)/(64*w**3) - 675675*I*t*exp(-I*t*w)/(16*w**5) - 230945*I*t**9*exp(-I*t*w)/(256*w) - 45045*I*t*exp(-I*t*w)/(128*w**3) - 45045*I*t**5*exp(-I*t*w)/(128*w) - 693*I*t*exp(-I*t*w)/(256*w) + 15015*I*t**3*exp(-I*t*w)/(256*w) + 88179*I*t**11*exp(-I*t*w)/(256*w) + 109395*I*t**7*exp(-I*t*w)/(128*w) + 225225*I*t**3*exp(-I*t*w)/(32*w**3) + 2078505*I*t**7*exp(-I*t*w)/(32*w**3) + 11486475*I*t**3*exp(-I*t*w)/(16*w**5) + 43648605*I*t**7*exp(-I*t*w)/(16*w**5) + 218243025*I*t**3*exp(-I*t*w)/(4*w**7) + 4583103525*I*t**3*exp(-I*t*w)/(2*w**9))*exp(I*t*w)]



simplifiedDirect=[
+ I
,
+ 1
+ I*w*t
,
- r3*I
+ r3*t*w

- r1/r2/f0*I*w2
+ r3/r0/f2*I*t2*w2
,
- r5/r0/f0

- r5/r0/f1*I*t*w

- r3/r2/f0*w2
+ r5/r0/f2*t2*w2

- r3/r2/f1*I*t*w3
+ r5/r0/f3*I*t3*w3
,
+ r7/r0/f0*I

- r7/r0/f1*t*w

+ r5/r2/f0*I*w2
- r7/r0/f2*I*t2*w2

- r5/r2/f1*t*w3
+ r7/r0/f3*t3*w3

+ r3/r4/f0*I*w4
- r5/r2/f2*I*t2*w4
+ r7/r0/f4*I*t4*w4
,
+ r9
+ r9*I*t*w
+ r7/b1*w2
- r9/b1*t2*w2
+ r7/b1*I*t*w3
- r9/r3/b1*I*t3*w3
+ r5/b3*w4
- r7/b2*t2*w4
+ r9/r3/b3*t4*w4
+ r5/b3*I*t*w5
- r7/r3/b2*I*t3*w5
+ r9/r5/b3*I*t5*w5
,
- r11*I
+ r11*t*w
- r9/b1*I*w2
+ r11/b1*I*t2*w2
+ r9/b1*t*w3
- r11/r3/b1*t3*w3
- r7/b3*I*w4
+ r9/b2*I*t2*w4
- r11/r3/b3*I*t4*w4
+ r7/b3*t*w5
- r9/r3/b2*t3*w5
+ r11/r5/b3*t5*w5
- r5/r3/b4*I*w6
+ r7/b4*I*t2*w6
- r9/r3/b4*I*t4*w6
+ r11/r5/r3/b4*I*t6*w6
,
- r13
- r11/b1*w2
- r9/b3*w4
- r7/r3/b4*w6
- r13/r3/b3*t4*w4
- r11/r3/b4*t4*w6
+ r9/b4*t2*w6
+ r13/r5/r3/b4*t6*w6
+ r11/b2*t2*w4
+ r13/b1*t2*(w2)
- r13*I*t*w
- r11/b1*I*t*w3
- r13/r5/b3*I*t5*w5
- r9/b3*I*t*w5
- r7/r3/b4*I*t*w7
+ r9/r3/b4*I*t3*w7
- r11/r5/b4*I*t5*w7
+ r13/r7/r3/b4*I*t7*w7
+ r11/r3/b2*I*t3*w5
+ r13/r3/b1*I*t3*w3
,
- r15*t*w
+ r15*I

+ r13/b1*I*w2
- r15/b1*I*t2*w2

- r13/b1*t*w3
+ r15/r3/b1*t3*w3

+ r11/b3*I*w4
- r13/b2*I*t2*w4
+ r15/r3/b3*I*t4*w4

- r11/b3*t*w5
+ r13/r3/b2*t3*w5
- r15/r5/b3*t5*w5

+ r9/r3/b4*I*w6
- r11/b4*I*t2*w6
+ r13/r3/b4*I*t4*w6
- r15/r5/r3/b4*I*t6*w6

- r9/r3/b4*t*w7
+ r11/r3/b4*t3*w7
- r13/r5/b4*t5*w7
+ r15/r7/r3/b4*t7*w7

+ r7/r3/b7*I*w8
- r9/r3/b5*I*t2*w8
+ r11/r3/b6*I*t4*w8
- r13/r5/r3/b5*I*t6*w8
+ r15/r7/r3/b7*I*t8*w8

,
+ r17

+ r17*I*t*w

+ r15/b1*w2
- r17/b1*t2*w2

+ r15/b1*I*t*w3
- r17/r3/b1*I*t3*w3

+ r13/b3*w4
- r15/b2*t2*w4
+ r17/r3/b3*t4*w4

+ r13/b3*I*t*w5
- r15/r3/b2*I*t3*w5
+ r17/r5/b3*I*t5*w5

+ r11/r3/b4*w6
- r13/b4*t2*w6
+ r15/r3/b4*t4*w6
- r17/r5/r3/b4*t6*w6

+ r11/r3/b4*I*t*w7
- r13/r3/b4*I*t3*w7
+ r15/r5/b4*I*t5*w7
- r17/r7/r3/b4*I*t7*w7

+ r9/r3/b7*w8
- r11/r3/b5*t2*w8
+ r13/r3/b6*t4*w8
- r15/r5/r3/b5*t6*w8
+ r17/r7/r3/b7*t8*w8

+ r9/r3/b7*I*t*w9
- r11/r3/r3/b5*I*t3*w9
+ r13/r5/b6*I*t5*w9
- r15/r7/r3/b5*I*t7*w9
+ r17/r9/r3/b7*I*t9*w9
,
- r19*I
+ r19*t*w
+ r19*I*t2*w2/b1
- r17*I*w2/b1
- r19/r3*t3*w3/b1
+ r17*t*w3/b1
- r15*I*w4/b3
+ r17*I*t2*w4/b2
- r19/r3*I*t4*w4/b3
+ r15*t*w5/b3
- r17/r3*t3*w5/b2
+ r19/r5*t5*w5/b3
- r13/r3*I*w6/b4
+ r15*I*t2*w6/b4
+ r19/r5/r3*I*t6*w6/b4
- r17/r3*I*t4*w6/b4

+ r13/r3/b4*t*w7
- r15/f3/b3*t3*w7
+ r17/f5/b1*t5*w7
- r19/f7*t7*w7

- r11/f8*r7*I*w8
+ r13/r3/b5*I*t2*w8
- r15/f4/b3*I*t4*w8
+ r17/f6/b1*I*t6*w8
- r19/f8*I*t8*w8

+ r11/r3*t*w9/b7
- r13/r3/r3/b5*t3*w9
+ r15/f5/b3*t5*w9
- r17/f7/b1*t7*w9
+ r19/f9*t9*w9

- r9/r5*I*w10/b8
+ r11/r3*I*t2*w10/b8
- r13/r3/r3*I*t4*w10/b7
+ r13/r3*I*t6*w10/b7
- r17/r7/r3*I*t8*w10/b8
+ r19/r0/f10*I*t10*w10
,
- r21
- r21/f1*I*t*w
- r19/f1/b1*w2
+ r21/f2*t2*w2
- r19/f1/b1*I*t*w3
+ r21/f3*I*t3*w3
- r17/f2/b2*w4
+ r19/f2/b1*t2*w4
- r21/f4*t4*w4
- r17/f1/b3*I*t*w5
+ r19/f3/b1*I*t3*w5
- r21/f5*I*t5*w5
- r15/f3/b3*w6
+ r17/f2/b3*t2*w6
- r19/f4/b1*t4*w6
+ r21/f6*t6*w6
- r15/f3/b3*I*t*w7
+ r17/r3/b4*I*t3*w7
- r19/r5/b4*I*t5*w7
+ r21/f7*I*t7*w7

- r13/f4/b4*w8
+ r15/f3/b4*t2*w8
- r17/f4/b3*t4*w8
+ r19/f6/b1*t6*w8
- r21/f8*t8*w8

- r13/r8/f1*I*t*w9
+ r15/r6/f3*I*t3*w9
- r17/r4/f5*I*t5*w9
+ r19/r2/f7*I*t7*w9
- r21/r0/f9*I*t9*w9

- r11/r10/f0*w10
+ r13/r8/f2*t2*w10
- r15/r6/f4*t4*w10
+ r17/r4/f6*t6*w10
- r19/r2/f8*t8*w10
+ r21/r0/f10*t10*w10

- r11/r10/f1*I*t*w11
+ r13/r8/f3*I*t3*w11
- r15/r6/f5*I*t5*w11
+ r17/r4/f7*I*t7*w11
- r19/r2/f9*I*t9*w11
+ r21/r0/f11*I*t11*w11
,
]

formulaDirect = [ I**(1-n) *
	sum([
		I**wi
		* w**wi 
		* sum([
			(-1)**( wi//2-i )
			* simplify(
				R(2*n -1 -wi+ti)
				/ R(wi-ti)
				/ factorial(ti)
			)
			* t**ti
			for i, ti in zip(xrange(wi//2+1), xrange(wi&1, wi+1, 2))
		])
		for wi in xrange((n+1))
	])
	for n in xrange(12)
]

# Test: should be [0,0,0...]
print [simplify(target_direct[i] - formulaDirect[i]) for i in xrange(8)]


full = [
	simplify(
	(2*i+1)**(1/2)/2 *( # TODO and still misses fft factor 1/2pi
		+ sh.subs(t,sympy.pi/2)*E.subs(t,sympy.pi/2) - sh.subs(t,0)
		+ di.subs(t,0)/w**(i+1) - di.subs(t,-1)*E.subs(t,-1)/w**(i+1)
	))
	for i, (sh, di) in enumerate(zip(simplifiedShadow, formulaDirect))
]

for i,f in enumerate(full): print i,":",f

import headDistortion
print "plotting..."
import numpy, os
import bmaudio
plot = bmaudio.SpectrumDisplay()
plot.inDb()
#plot.setLogFrequency()
plot.ylim(-35,5)
ws=numpy.arange(0,22050*2*numpy.pi,spectralRange*2*numpy.pi/spectrumBins*20)
for i,f in enumerate(full[:7]) :
	print i
#	plot.addSpectrumData(
#		headDistortion.sphericalHead3dSHSpectrum(i, spectralRange, spectrumBins/20, to*c, c),
#		22050, "H%i,0,+"%i)
	plot.addSpectrumData(
		numpy.array([ complex(f.subs(w,wi*to)) for wi in ws ]),
		22050, "H%i,0,+"%i)
plot.hardcopy('figures/'+os.path.basename(os.path.splitext(__file__)[0])+".pdf", "pdf")







