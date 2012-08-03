#!/usr/bin/env python

from spectral_components3d_common import *
from spectral_components3d_Pli import Pli

"""
The problem findng the spectral domain expression of Hl,0,+:
sqrt(2*l+1)/2/to * ( # TODO: maybe some 1/2pi is needed for fourier but this makes H0(0) be 0dB
	integrate(Pl(t/to) * exp(-I*w*t), (t,-to,0)) +
	integrate(cos(t/to) * Pl(sin(t/to)) * exp(-I*w*t), (t,0,to*pi/2)) 
	)
	Substitution:
		t = to * z
		dt = to * dz
		z0=-1 z1=0
		z0=0, z1= pi/2
sqrt(2*l+1)/2 * (
	+ integrate(Pl(z)*exp(-1j*w*to*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(sin(z))*exp(-1j*w*to*z), (z,0,pi/2) # shadowed part
	)
	Substitution:
		y = -1j*w*to
sqrt(2*l+1)/2 * (
	+ integrate(Pl(z)*exp(y*z), (z,-1, 0))             # direct part
	+ integrate(cos(z)*Pl(sin(z))*exp(y*z), (z,0,pi/2) # shadowed part
	)
"""

"""
################################################################
# Direct part integral
################################################################
So now we focus on the direct part
	integrate(Pl(z)*exp(y*z), (z,-1, 0))
"""


if __name__ == "__main__" :

	maz.check("DirectIndefiniteIntegralOnZY", "start z=t/to, y=-I*w*to, %i", [
		integrate(Pl(l,z) * exp(y*z), z)
	for l in xrange(10)])
	maz.check("DirectIndefiniteIntegralOnZY", "apply recursive integration by parts, %i", [
		sum([
			exp(y*z) * (-1)**N(i) * Pli(l,i,z) * y**N(-i-1)
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "move exp outside, %i", [
		exp(y*z) * 
		sum([
			(-1)**N(i) * Pli(l,i,z) * y**N(-i-1)
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "inline Pli, %i", [
		exp(y*z) * 
		sum([
			(-1)**N(i) 
			* y**N(-i-1)
			* sum([
				(-1)**N(k//2)
				/ N(2)**N(l)
				* f(2*l-k)
				/ f(l-k//2)
				/ f(l-k-i)
				/ f(k//2)
				* z**N(l-k-i)
			for k in xrange(0,l-i+1,2) ])
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "all factors inside, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k//2+i)
			/ N(2)**N(l)
			* f(2*l-k)
			/ f(l-k//2)
			/ f(l-k-i)
			/ f(k//2)
			* y**N(-i-1)
			* z**N(l-k-i)
		for k in xrange(0,l-i+1,2) ])
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "double k, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
		for k in xrange((l-i)//2+1) ])
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "sum interval into conditions, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
			if True
			and 2*k <= l-i
			and 2*k >= 0
			and i <= l+1
			and i >=0
			else 0
		for k in xrange((l-i)//2+1) ])
		for i in xrange(l+2)])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "inverting vars, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
			if True
			and i >=0
			and i <= ( l-2*k if l-2*k <= l+1 else l+1)
			and 2*k >= 0
			else 0
		for i in xrange(l+2)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "split summatory vars, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
			if True
			and i >=0
			and i <= l-2*k
			and -2*k <= 1
			and 2*k >= 0
			else 0
		for i in xrange(l+2)])
		+
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
			if True
			and i >=0
			and i <= l+1
			and 2*k < 0 # note that those are complementary
			and 2*k >= 0
			else 0
		for i in xrange(l+2)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "second summatory is empty, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
			if True
			and i >=0
			and i <= l-2*k
			and -2*k <= 1
			and 2*k >= 0
			else 0
		for i in xrange(l+2)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "second summatory is empty and first i bounds changed, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-2*k-i)
			/ f(k)
			* y**N(-i-1)
			* z**N(l-2*k-i)
		for i in xrange(l-2*k+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectIndefiniteIntegralOnZY", "i'=i+2*k, %i", [
		exp(y*z) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* z**N(l-i)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "z in -1,0 interval, %i", [
		eq.subs(z,0) - eq.subs(z,-1)
	for eq in maz.recall("DirectIndefiniteIntegralOnZY") ])

	maz.check("DirectDefiniteIntegralOnY", "inline functions and subsitute, %i", [
		exp(y*N(0)) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* N(0)**N(l-i)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
		-
		exp(y*N(-1)) * 
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* N(-1)**N(l-i)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])
	maz.check("DirectDefiniteIntegralOnY", "move exponential inside, %i", [
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* N(0)**N(l-i)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
		-
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* N(-1)**N(l-i)
			/ exp(y)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "first term is non zero just when i=l, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(-l+2*k-1)
		for k in xrange(l//2+1) ])
		-
		sum([
		sum([
			(-1)**N(k+i)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			* N(-1)**N(l-i)
			/ exp(y)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "first term is non zero just when i=l, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(-l+2*k-1)
		for k in xrange(l//2+1) ])
		-
		sum([
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(l-i)
			/ f(k)
			* y**N(-i+2*k-1)
			/ exp(y)
		for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "joining k summatories, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(-l+2*k-1)
			-
			sum([
				(-1)**N(k+l)
				/ N(2)**N(l)
				* f(2*l-2*k)
				/ f(l-k)
				/ f(l-i)
				/ f(k)
				* y**N(-i+2*k-1)
				/ exp(y)
			for i in xrange(2*k,l+1)])
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "extract common factors, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(-l+2*k-1)
			* (
			1
			-
			sum([
				y**N(-i+2*k-1)
				/ y**N(-l+2*k-1)
				/ f(l-i)
			for i in xrange(2*k,l+1)])
			/ exp(y)
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "simplify y power, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(2*k-l-1)
			* (
			1
			-
			sum([
				y**N(l-i)
				/ f(l-i)
			for i in xrange(2*k,l+1)])
			/ exp(y)
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "i'=l-i, %i", [
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(2*k-l-1)
			* (
			1
			-
			sum([
				y**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
			/ exp(y)
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegralOnY", "i'=l-i, %i", [
		exp(-y) *
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* y**N(2*k-l-1)
			* (
			exp(y) -
			sum([
				y**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	for formula in maz.recall("DirectDefiniteIntegralOnY") : print simplify(formula); print


	maz.check("DirectDefiniteIntegral", "y=-I*w*to, %i", [
		eq.subs(y,-I*w*to)
	for eq in maz.recall("DirectDefiniteIntegralOnY") ])

	maz.check("DirectDefiniteIntegral", "inlining, %i", [
		exp(I*w*to) *
		sum([
			(-1)**N(k+l)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (-I*w*to)**N(2*k-l-1)
			* (
			exp(-I*w*to) -
			sum([
				(-I*w*to)**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "oscillators out, %i", [
		exp(I*w*to) *
		sum([
			I**N(-l+1)
			/ N(2)**N(l)
			* f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (w*to)**N(2*k-l-1)
			* (
			exp(-I*w*to) -
			sum([
				I**N(-i)
				* (w*to)**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "common factors out, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* sum([
			f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (w*to)**N(2*k-l-1)
			* (
			1 - exp(I*w*to)*
			sum([
				I**N(-i)
				* (w*to)**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
			)
		for k in xrange(l//2+1) ])
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "split terms, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (w*to)**N(2*k-l-1)
		for k in xrange(l//2+1) ])
		- exp(I*w*to) * sum([
			f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (w*to)**N(2*k-l-1)
			*
			sum([
				I**N(-i)
				* (w*to)**N(i)
				/ f(i)
			for i in xrange(l-2*k+1)])
		for k in xrange(l//2+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "introduce new vars, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
		sum([
			f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			* (w*to)**N(2*k-l-1)
			if True
			and m==2*k-l-1
			else 0
		for k in xrange(l//2+1) ])
		for m in xrange(-l-1,l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(2*l-2*k)
			/ f(l-k)
			/ f(k)
			*
			sum([
				I**N(-i)
				* (w*to)**N(2*k+i-l-1)
				/ f(i)
			if True
			and m==2*k+i-l-1
			else 0
			for i in xrange(l-2*k+1)])
		for k in xrange(l//2+1) ])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "k'=2*k, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			* (w*to)**N(k-l-1)
			if True
			and m==k-l-1
			else 0
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			*
			sum([
				I**N(-i)
				* (w*to)**N(k+i-l-1)
				/ f(i)
			if True
			and m==k+i-l-1
			else 0
			for i in xrange(l-k+1)])
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "k bounds to conditions, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			* (w*to)**N(k-l-1)
			if True
			and m==k-l-1
			and k<=l
			and k >= 0
			else 0
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			*
			sum([
				I**N(-i)
				* (w*to)**N(k+i-l-1)
				/ f(i)
			if True
			and m==k+i-l-1
			and k <= l
			and k >= 0
			else 0
			for i in xrange(l-k+1)])
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "substitute k=m+l+1, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l-m-1)
			/ f((l-m-1)//2)
			/ f((l+m+1)//2)
			* (w*to)**N(m)
			if True
			and m <= -1
			and m >= -l-1
			and (m+l+1)%2 == 0
			else 0
		for m in xrange(-l-1,0) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			*
			sum([
				I**N(-i)
				* (w*to)**N(k+i-l-1)
				/ f(i)
			if True
			and m==k+i-l-1
			and k <= l
			and k >= 0
			else 0
			for i in xrange(l-k+1)])
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "substitute m'=-m-1, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			if True
			and (-m+l)%2 == 0
			else 0
		for m in xrange(0,l+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			*
			sum([
				I**N(-i)
				* (w*to)**N(k+i-l-1)
				/ f(i)
			if True
			and m==k+i-l-1
			and k <= l
			and k >= 0
			else 0
			for i in xrange(l-k+1)])
		for k in xrange(0,l+1,2) ])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "swap i and k by relaxing i and adding a condition, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			if True
			and (-m+l)%2 == 0
			else 0
		for m in xrange(0,l+1) ])
		- exp(I*w*to) *
		sum([
		sum([
		sum([
			f(2*l-k)
			/ f(l-k//2)
			/ f(k//2)
			* I**N(-i)
			* (w*to)**N(k+i-l-1)
			/ f(i)
			if True
			and m==k+i-l-1
			and k <= l
			and k >= 0
			and i <= l-k
			else 0
		for k in xrange(0,l+1,2) ])
		for i in xrange(l+1)])
		for m in xrange(-l-1,l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "k=m+l+1-i, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			if True
			and (-m+l)%2 == 0
			else 0
		for m in xrange(0,l+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(-m+l-1+i)
			/ f((-m+l-1+i)//2)
			/ f((m+l+1-i)//2)
			* I**N(-i)
			/ f(i)
			if True
			and m+1 <= i
			and l+m+1 >= i
			and (m+l+1-i)%2 == 0
			else 0
		for i in xrange(l+1)])
		* (w*to)**N(m)
		for m in xrange(-l-1,0) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "m'=-m-1, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			if True
			and (l-m)%2 == 0
			else 0
		for m in xrange(l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(l+m+i)
			* I**N(-i)
			/ f((l+m+i)//2)
			/ f((l-m-i)//2)
			/ f(i)
			if True
			and (l-m-i)%2 == 0
			else 0
		for i in xrange(l-m+1)])
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "i'=i+m, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			if True
			and (l-m)%2 == 0
			else 0
		for m in xrange(l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(l+i)
			* I**N(-i+m)
			/ f((l+i)//2)
			/ f((l-i)//2)
			/ f(i-m)
			if True
			and (l-i)%2 == 0
			else 0
		for i in xrange(m,l+1)])
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "join summatories, %i", [
		I**N(-l+1)
		/ N(2)**N(l)
		* (
		sum([
			f(l+m)
			/ f((l-m)//2)
			/ f((l+m)//2)
			/ (w*to)**N(m+1)
			* (N(1)+N(-1)**(l-m))
			/ N(2)
		for m in xrange(l+1+1) ])
		- exp(I*w*to) *
		sum([
		sum([
			f(l+i)
			* I**N(-i+m)
			/ f((l+i)//2)
			/ f((l-i)//2)
			/ f(i-m)
			* (N(1)+N(-1)**(l-i))
			/ N(2)
		for i in xrange(m,l+1)])
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "join summatories, %i", [
		sum([
		(
			(N(1)+N(-1)**(l-m))
			* I**N(-m)
			* f(l+m)
			/ f((l+m)//2)
			/ f((l-m)//2)
			- exp(I*w*to) *
			sum([
				(N(1)+N(-1)**(l-i))
				* f(l+i)
				* I**N(-i)
				/ f((l+i)//2)
				/ f((l-i)//2)
				/ f(i-m)
			for i in xrange(m,l+1)])
		)
		* I**N(m)
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		* I**N(-l+1)
		/ N(2)**N(l+1)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "Use double factorials, %i", [
		sum([
		(
			(N(1)+N(-1)**(l-m))
			* I**N(-m)
			* R(l+m-1)
			* R(l-m-1)
			/ f(l-m)
			- exp(I*w*to) *
			sum([
				(N(1)+N(-1)**(l-i))
				* I**N(-i)
				* R(l+i-1)
				* R(l-i-1)
				/ f(l-i)
				/ f(i-m)
			for i in xrange(m,l+1)])
		)
		* I**N(m)
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		* I**N(-l+1)
		/ N(2)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "-i'=i-l, %i", [
		sum([
		(
			(N(1)+N(-1)**(l-m))
			* I**N(-m)
			* R(l+m-1)
			* R(l-m-1)
			/ f(l-m)
			- exp(I*w*to) *
			sum([
				(N(1)+N(-1)**(i))
				* I**N(i-l)
				* R(2*l-i-1)
				* R(i-1)
				/ f(i)
				/ f(+l-i-m)
			for i in xrange(0,l-m+1)])
		)
		* I**N(m)
		/ (w*to)**N(m+1)
		for m in xrange(l+1+1) ])
		* I**N(-l+1)
		/ N(2)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "-m'=m-l, %i", [
		sum([
		(
			(N(1)+N(-1)**(m))
			* I**N(m-l)
			* R(2*l-m-1)
			* R(m-1)
			/ f(m)
			- exp(I*w*to) *
			sum([
				(N(1)+N(-1)**(i))
				* I**N(i-l)
				* R(2*l-i-1)
				* R(i-1)
				/ f(i)
				/ f(m-i)
			for i in xrange(0,m+1)])
		)
		* I**N(l-m)
		/ (w*to)**N(l-m+1)
		for m in xrange(0,l+1) ])
		* I**N(-l+1)
		/ N(2)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "R(k-1)/f(k)=R(k+1)/f(k+1)=1/R(k), %i", [
		sum([
		(
			(N(1)+N(-1)**(m))
			* I**N(m-l)
			* R(2*l-m-1)
			/ R(m)
			- exp(I*w*to) *
			sum([
				(N(1)+N(-1)**(i))
				* I**N(i-l)
				* R(2*l-i-1)
				/ R(i)
				/ f(m-i)
			for i in xrange(0,m+1)])
		)
		* I**N(l-m)
		/ (w*to)**N(l-m+1)
		for m in xrange(0,l+1) ])
		* I**N(-l+1)
		/ N(2)
	for l in xrange(12)])

	maz.check("DirectDefiniteIntegral", "common factors on I, %i", [
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
		* I**N(-m)
		/ (w*to)**N(l-m+1)
		for m in xrange(l+1) ])
		* I**N(-l+1)
	for l in xrange(12)])


def DirectDefiniteIntegral(l,w) :
	return (
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

if __name__ == "__main__" :

	maz.check("DirectDefiniteIntegral", "Function, %i", [
		DirectDefiniteIntegral(l,w)
	for l in xrange(12)])


	for formula in maz.recall("DirectDefiniteIntegral") : print formula; print




