#!/usr/bin/env python
from __future__ import division
from spectral_components3d_common import *

# Pli is the ith derivative of the Legendre polynomial of degre l

"""
# Previous problem: ith derivative of the legendre polynomial of l degree Pil(l,i,x)

	BN(i,j) = binomial # = f(i)/f(j)/f(i-j)
Pl(x)     = 2**(-l) * sum( [ (-1)**(k) * BN(l,k) * BN(2*(l-k), l)           * x**(l-2*k) for k in xrange(l//2+1) ] )
d1(Pl(x)) = 2**(-l) * sum( [ (-1)**(k) * BN(l,k) * BN(2*(l-k), l) * f(l-2*k) * x**(l-2*k-1) for k in xrange((l-1)//2+1) ] )
di(Pl(x)) = 2**(-l) * sum( [ (-1)**(k) * BN(l,k) * BN(2*(l-k), l) * f(l-2*k) / f(l-2*k-i) * x**(l-2*k-i) for k in xrange((l-i)//2+1) ] )
				2**(-l) * BN(l,k) * BN(2*(l-k), l) * f(l-2*k) / f(l-2*k-i) =
					= 2**(-l) * f(l) / f(k) / f(l-k) * f(2*(l-k)) / f(l) / f(l-2*k) * f(l-2*k) / f(l-2*k-i)
					= 2**(-l) / f(k) / f(l-k) * f(2*(l-k)) / f(l-2*k) * f(l-2*k) / f(l-2*k-i)
					= 2**(-l) / f(k) / f(l-k) * f(2*(l-k)) / f(l-2*k-i)
					= 2**(-l) * f(2*(l-k)) / f(l-k) / f(l-2*k-i) / f(k)
						f(2n)/f(n) = 2**n * R(2n-1)
					= 2**((l-k)-l) * R(2*(l-k)-1) / f(l-2*k-i) / f(k)
					= 2**(-k) * R(2*l-2*k-1) / f(l-2*k-i) / f(k)
				
	= sum( [ (-1)**(k) * 2**(-k) * R(2*l-2*k-1) / f(l-2*k-i) / f(k) * x**(l-2*k-i) for k in xrange((l-i)//2+1) ] )
	= sum( [ (-2)**(-k) * R(2*l-2*k-1) / f(l-2*k-i) / f(k) * x**(l-2*k-i) for k in xrange((l-i)//2+1) ] )
	= sum( [ (-2)**(-k) * R(2*(l-k)-1) / f(l-2*k-i) / f(k) * x**(l-2*k-i) for k in xrange((l-i)//2+1) ] )
	= sum( [
		(I**(l-i-k)+I**(-(l-i-k)))/2
		/ 2**((l-i-k)//2)
		* R(l+i+k-1)
		/ f(k)
		/ f((l-i-k)//2)
		* x**k
		for k in xrange(l-i+1) ] )

Check: For some orders compare the sympy result with the formula we got.
"""

def Pli_bysympy(l,i,x) : return sympy.diff(Pl(l,x),x,i)

def Pli(l,i,x) :
	"""ith derivative of the l order legendre polynomial"""
	# doubling k
	return sum([
		(-1)**N(k//2)
		/ N(2)**N(l)
		* f(2*l-k)
		/ f(l-k//2)
		/ f(l-k-i)
		/ f(k//2)
		* x**N(l-k-i)
		for k in xrange(0,l-i+1,2) ] )
	# R(2*(l-k)-1) = f(2*(l-k)) / 2**(l-k) / f(l-k)
	return sum([
		(-1)**N(k)
		/ N(2)**N(l)
		* f(2*l-2*k)
		/ f(l-k)
		/ f(l-2*k-i)
		/ f(k)
		* x**N(l-2*k-i)
		for k in xrange(0,(l-i+2)//2) ] )
	# k' = 2k
	return sum([
		(-1)**k
		/ 2**(k)
		* R(2*l-2*k-1)
		/ f(-2*k-i+l)
		/ f(k)
		* x**(-2*k-i+l)
		for k in xrange(0,(l-i+2)//2) ] )
	# k' = -k
	return sum([
		(I**(-k)+I**(k))/2
		/ 2**(k//2)
		* R(2*l-k-1)
		/ f(-k-i+l)
		/ f(k//2)
		* x**(-k-i+l)
		for k in xrange(0,l-i+1) ] )
	# simplified
	return sum([
		(I**(-k)+I**(k))/2
		/ 2**(-k//2)
		* R(2*l+k-1)
		/ f(k-i+l)
		/ f(-k//2)
		* x**(k-i+l)
		for k in xrange(0,-1+i-l,-2) ] )
	# k' = k-l+i
	return sum([
		(I**(l-i-(k-i+l))+I**(-(l-i-(k-i+l))))/2
		/ 2**((l-i-(k-i+l))//2)
		* R(l+i+(k-i+l)-1)
		/ f((k-i+l))
		/ f((l-i-(k-i+l))//2)
		* x**(k-i+l)
		for k in xrange(0,-1+i-l,-2) ] )
	return sum([
		(I**(l-i-k)+I**(-(l-i-k)))/2
		/ 2**((l-i-k)//2)
		* R(l+i+k-1)
		/ f(k)
		/ f((l-i-k)//2)
		* x**k
		for k in xrange(l-i,-1,-2) ] )
	return sum([
		(I**(l-i-k)+I**(-(l-i-k)))/2
		/ 2**((l-i-k)//2)
		* R(l+i+k-1)
		/ f(k)
		/ f((l-i-k)//2)
		* x**k
		for k in xrange(l-i+1) ] )

if __name__ == "__main__" :

	cases = reduce(operator.add,
		[[ (l,i)
			for i in xrange(l+3) ]
			for l in xrange(7) ],
		[])
	maz.cases("Pli", cases)

	maz.check("Pli", "Computed by sympy one by one, %s", [
		Pli_bysympy(l,i,x)
	for l,i in cases])

	maz.check("Pli", "General formula, %s", [
		Pli(l,i,x)
	for l,i in cases])





