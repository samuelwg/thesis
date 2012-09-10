#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Copyright 2012 David García Garzón

This file is part of spherelab

spherelab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

spherelab is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from PySide import QtCore, QtGui, QtOpenGL
import sys
import numpy as np
import time
import math

order = 5
shShape = (order+1,order+1)
shSize = np.prod(shShape)


def shIndexes(order) :
	"""Returns an enumeration of SH order and degrees"""
	return [(l,m) for l in xrange(order+1) for m in xrange(-l,l+1) ]

def shi(l,m) :
	"""Returns the mapping of SH order (l) and degree (m) to a square matrix index"""
	return (l-m,l) if m>0 else (l,l+m)

def shi_reverse(i,j) :
	l = max(i,j)
	absm = l-min(i,j)
	return (l,absm) if i<j else (l,-absm)

def ead2xyz(e,a,d) :
	"""
		Converts polars to cartesian with those conventions:
		x front, y left, z up
		0 elevation and azimuth front
		positive elevation up, positive azimuth left
	"""
	ra, re = np.radians(a), np.radians(e)
	sa, se = np.sin(ra), np.sin(re)
	ca, ce = np.cos(ra), np.cos(re)
	x,y,z = d*ce*ca, d*ce*sa, d*se
	return x,y,z


########################################
# Normalization factors,
# They are to be applied to the sn3d to get a different normalization.


# factors to be applied to a sn3d sh to get n3d normalization.
# n3d normalization is that so that they have unit power
n3d = np.zeros(shShape)
for l in xrange(order+1) :
	for m in xrange(-l, l+1) :
		n3d[shi(l,m)] = math.sqrt(2*l+1)

# real orthonormalized with an an extra factor of 1/sqrt(4 pi)
# Ensure that the autoconvolution is one instead of 1/4pi
on3d = n3d/np.sqrt(4*np.pi)

# factors to be applied to sn3d to get maxn normalization.
# maxn normalization is that so that absolute maximum value
# for each sh is 1.
maxn=np.zeros(shShape)
maxn[shi(0, 0)] = 1.

maxn[shi(1,+1)] = 1.
maxn[shi(1, 0)] = 1.
maxn[shi(1,-1)] = 1.

maxn[shi(2,+2)] = math.sqrt(4./3)
maxn[shi(2,+1)] = math.sqrt(4./3)
maxn[shi(2, 0)] = 1.
maxn[shi(2,-1)] = math.sqrt(4./3)
maxn[shi(2,-2)] = math.sqrt(4./3)

maxn[shi(3,+3)] = math.sqrt(8./5)
maxn[shi(3,+2)] = math.sqrt(9./5)
maxn[shi(3,+1)] = math.sqrt(45./32)
maxn[shi(3, 0)] = 1.
maxn[shi(3,-1)] = math.sqrt(45./32)
maxn[shi(3,-2)] = math.sqrt(9./5)
maxn[shi(3,-3)] = math.sqrt(8./5)

# TODO: maxn for 4th order and beyond
maxn[shi(4,+4)] = 1.
maxn[shi(4,+3)] = 1.
maxn[shi(4,+2)] = 1.
maxn[shi(4,+1)] = 1.
maxn[shi(4, 0)] = 1.
maxn[shi(4,-1)] = 1.
maxn[shi(4,-2)] = 1.
maxn[shi(4,-3)] = 1.
maxn[shi(4,-4)] = 1.

maxn[shi(5,+5)] = 1.
maxn[shi(5,+4)] = 1.
maxn[shi(5,+3)] = 1.
maxn[shi(5,+2)] = 1.
maxn[shi(5,+1)] = 1.
maxn[shi(5, 0)] = 1.
maxn[shi(5,-1)] = 1.
maxn[shi(5,-2)] = 1.
maxn[shi(5,-3)] = 1.
maxn[shi(5,-4)] = 1.
maxn[shi(5,-5)] = 1.

# factors to be applied to sn3d sh to get fuma normalization.
# fuma is like maxn but applying a sqrt(1./2) factor to the 0,0 channel
# in order to be compatible with B-Format standard for the first order.
fuma=maxn.copy()
fuma[shi(0, 0)] = math.sqrt(1./2)


class SemiNormalizedSH(object) :

	def __init__(self, compiled=False) :
		from sympy import S, symbols, sqrt, Lambda, lambdify, Matrix
		from sympy.utilities.autowrap import autowrap
		x,y,z = symbols("x y z")

		self.order = order
		self.sn3d = np.array([[None]*(order+1)]*(order+1))

		self.sn3d[shi(0, 0)] = S(1)

		self.sn3d[shi(1,+1)] = x
		self.sn3d[shi(1, 0)] = z
		self.sn3d[shi(1,-1)] = y

		self.sn3d[shi(2,-2)] = (x*y)*2       * sqrt(S(3)/4)
		self.sn3d[shi(2,-1)] = (y*z)*2       * sqrt(S(3)/4)
		self.sn3d[shi(2, 0)] = (3*z*z -1)/2
		self.sn3d[shi(2,+1)] = (x*z)*2       * sqrt(S(3)/4)
		self.sn3d[shi(2,+2)] = (x*x - y*y)   * sqrt(S(3)/4)

		self.sn3d[shi(3,+3)] = x*(x*x-3*y*y) * sqrt(S(5)/8)
		self.sn3d[shi(3,+2)] = z*(x*x-y*y)   * sqrt(S(15)/4)
		self.sn3d[shi(3,+1)] = x*(5*z*z -1)  * sqrt(S(3)/8)
		self.sn3d[shi(3, 0)] = z*(5*z*z -3)  * sqrt(S(1)/4)
		self.sn3d[shi(3,-1)] = y*(5*z*z -1)  * sqrt(S(3)/8)
		self.sn3d[shi(3,-2)] = z*x*y*2       * sqrt(S(15)/4)
		self.sn3d[shi(3,-3)] = y*(3*x*x-y*y) * sqrt(S(5)/8)

		self.sn3d[shi(4,+4)] = (x**4 -6*x*x*y*y +y**4)  * sqrt(S(35)/64)
		self.sn3d[shi(4,+3)] = z*x*(x*x -3*y*y)         * sqrt(S(35)/8)
		self.sn3d[shi(4,+2)] = (x*x-y*y)*(7*z*z -1)     * sqrt(S(5)/16)
		self.sn3d[shi(4,+1)] = z*x*(7*z*z -3)           * sqrt(S(5)/8)
		self.sn3d[shi(4, 0)] = (3+z*z*(-30+z*z*35))     * sqrt(S(1)/64)
		self.sn3d[shi(4,-1)] = z*y*(7*z*z -3)           * sqrt(S(5)/8)
		self.sn3d[shi(4,-2)] = 2*x*y*(7*z*z -1)         * sqrt(S(5)/16)
		self.sn3d[shi(4,-3)] = z*y*(3*x*x -y*y)         * sqrt(S(35)/8)
		self.sn3d[shi(4,-4)] = 4*x*y*(x*x -y*y)         * sqrt(S(35)/64)

		self.sn3d[shi(5,+5)] = x*(x**4 -10*x*x*y*y +5*y**4) * sqrt(S(9)*7/64/2)
		self.sn3d[shi(5,+4)] = z*(x**4 -6*x*x*y*y +y**4)    * sqrt(S(35)*9/64)
		self.sn3d[shi(5,+3)] = x*(x*x-3*y*y)*(9*z*z-1)      * sqrt(S(5)*7/128)
		self.sn3d[shi(5,+2)] = 3*z*(x*x-y*y)*(3*z*z-1)      * sqrt(S(7)*5/3/16)
		self.sn3d[shi(5,+1)] = x*(21*z**4 -14*z**2 + 1)     * sqrt(S(3)*5/64)
		self.sn3d[shi(5, 0)] = z*(63*z**4 -70*z**2 + 15)    * sqrt(S(1)/64)
		self.sn3d[shi(5,-1)] = y*(21*z**4 -14*z**2 + 1)     * sqrt(S(3)*5/64)
		self.sn3d[shi(5,-2)] = 2*x*y*z*(3*z*z -1)           * sqrt(S(105)/16)
		self.sn3d[shi(5,-3)] = y*(9*z*z-1)*(3*x*x-y*y)      * sqrt(S(35)/8/16)
		self.sn3d[shi(5,-4)] = 4*x*y*z*(x*x -y*y)           * sqrt(S(9)*35./64)
		self.sn3d[shi(5,-5)] = y*(5*x**4 -7*x*x*y*y +y**4)  * sqrt(S(9)/10)
		self.sn3d[shi(5,-5)] = y*( 5*x**4 -10*x*x*y*y +y**4 )  * sqrt(S(9)*7/64/2)

		if compiled :
			# TODO: Get a single function that returns the whole matrix
			# Compiling a library for each cell takes soooo long
			# Hints for future work:
			# http://ojensen.wordpress.com/2010/08/10/fast-ufunc-ish-hydrogen-solutions/
			# https://groups.google.com/forum/?fromgroups=#!topic/sympy/KXoO-BvmD_w

			for i in xrange(self.order+1) :
				for j in xrange(self.order+1) :
					self.sn3d[i,j] = autowrap(
						self.sn3d[i,j],
						language='c',
						backend='cython',
						args=(x,y,z))
		else :
			for i in xrange(self.order+1) :
				for j in xrange(self.order+1) :
#					self.sn3d[i,j] = Lambda(
#						(x,y,z),
#						self.sn3d[i,j],
#						)
					self.sn3d[i,j] = lambdify(
						(x,y,z),
						self.sn3d[i,j],
						)

	def evalEA(self, e, a, target=None) :
		x,y,z = ead2xyz(e, a, 1)
		return self.evalXYZ(x,y,z, target)

	def evalXYZ(self, x,y,z, target=None) :
		if target is None : target = np.zeros(self.sn3d.shape)
		assert target.shape == self.sn3d.shape
		for i in xrange(self.order+1) :
			for j in xrange(self.order+1) :
				target[i,j] = self.sn3d[i,j](x,y,z)
		return target

print "Precomputing..."
from time import time
t = time()
obj = SemiNormalizedSH()
print "done, elapsed: %.4fs"%(time()-t)

def semiNormalizedSH_lambdified(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
	return obj.evalEA(e,a, target)

# TODO: Test compiled version


def semiNormalizedSH_genericNumpy(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation.
	Version based on the generic formula, using sympy just for legendre.
	"""
	x,y,z = ead2xyz(e, a, 1)

	sn3d = np.zeros((order+1,order+1)) if target is None else target
	from sympy import factorial, sqrt, S
	import sympy as sp
	ra = math.radians(a)

	for l in xrange(order+1) :
		for m in xrange(-l,l+1) :
			absm = abs(m)
			factor = math.sqrt(2) if m else 1.
			sn3d[shi(l,m)] = (factor
				/ math.sqrt(
					sp.factorial(l+absm) / sp.factorial(l-absm) # clearest for the next one
#					np.prod(xrange(l-absm+1, l+absm+1))
					)
				* ( sp.assoc_legendre(l,absm,z) * S(-1)**m ).evalf()
				* ( np.cos(m*ra) if m>=0 else np.sin(absm*ra) )
				)
	return sn3d

def semiNormalizedSH_sympyGeneratedExpressions(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	Version based on sympy expansion of the generic formula for each l,m.
	being a function of z (sin(elevation) and azimuth
	"""
	x,y,z = ead2xyz(e, a, 1)

	sn3d = np.zeros((order+1,order+1)) if target is None else target
	ra = math.radians(a)

	# sympy generated from generic formula
	sn3d[shi(0, 0)] = 1
	sn3d[shi(1,-1)] = (-z**2 + 1)**(1./2)*np.sin(ra)
	sn3d[shi(1, 0)] = z
	sn3d[shi(1,+1)] = (-z**2 + 1)**(1./2)*np.cos(ra)
	sn3d[shi(2,-2)] = 3**(1./2)*(-z**2 + 1)*np.sin(2*ra)/2
	sn3d[shi(2,-1)] = 3**(1./2)*z*(-z**2 + 1)**(1./2)*np.sin(ra)
	sn3d[shi(2, 0)] = 3*z**2/2 - 1./2
	sn3d[shi(2,+1)] = 3**(1./2)*z*(-z**2 + 1)**(1./2)*np.cos(ra)
	sn3d[shi(2,+2)] = 3**(1./2)*(-z**2 + 1)*np.cos(2*ra)/2
	sn3d[shi(3,-3)] = 10**(1./2)*(-z**2 + 1)**(3./2)*np.sin(3*ra)/4
	sn3d[shi(3,-2)] = 15**(1./2)*z*(-z**2 + 1)*np.sin(2*ra)/2
	sn3d[shi(3,-1)] = 6**(1./2)*(-z**2 + 1)**(1./2)*(5*z**2 - 1)*np.sin(ra)/4
	sn3d[shi(3, 0)] = z*(5*z**2 - 3)/2
	sn3d[shi(3,+1)] = 6**(1./2)*(-z**2 + 1)**(1./2)*(5*z**2 - 1)*np.cos(ra)/4
	sn3d[shi(3,+2)] = 15**(1./2)*z*(-z**2 + 1)*np.cos(2*ra)/2
	sn3d[shi(3,+3)] = 10**(1./2)*(-z**2 + 1)**(3./2)*np.cos(3*ra)/4
	sn3d[shi(4,-4)] = 35**(1./2)*(z**4 - 2*z**2 + 1)*np.sin(4*ra)/8
	sn3d[shi(4,-3)] = 70**(1./2)*z*(-z**2 + 1)**(3./2)*np.sin(3*ra)/4
	sn3d[shi(4,-2)] = 5**(1./2)*(-7*z**4 + 8*z**2 - 1)*np.sin(2*ra)/4
	sn3d[shi(4,-1)] = 10**(1./2)*z*(-z**2 + 1)**(1./2)*(7*z**2 - 3)*np.sin(ra)/4
	sn3d[shi(4, 0)] = 35*z**4/8 - 15*z**2/4 + 3./8
	sn3d[shi(4,+1)] = 10**(1./2)*z*(-z**2 + 1)**(1./2)*(7*z**2 - 3)*np.cos(ra)/4
	sn3d[shi(4,+2)] = 5**(1./2)*(-7*z**4 + 8*z**2 - 1)*np.cos(2*ra)/4
	sn3d[shi(4,+3)] = 70**(1./2)*z*(-z**2 + 1)**(3./2)*np.cos(3*ra)/4
	sn3d[shi(4,+4)] = 35**(1./2)*(z**4 - 2*z**2 + 1)*np.cos(4*ra)/8
	sn3d[shi(5,-5)] = 3*14**(1./2)*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.sin(5*ra)/16
	sn3d[shi(5,-4)] = 3*35**(1./2)*z*(z**4 - 2*z**2 + 1)*np.sin(4*ra)/8
	sn3d[shi(5,-3)] = 70**(1./2)*(-z**2 + 1)**(1./2)*(-9*z**4 + 10*z**2 - 1)*np.sin(3*ra)/16
	sn3d[shi(5,-2)] = 105**(1./2)*z*(-3*z**4 + 4*z**2 - 1)*np.sin(2*ra)/4
	sn3d[shi(5,-1)] = 15**(1./2)*(-z**2 + 1)**(1./2)*(21*z**4 - 14*z**2 + 1)*np.sin(ra)/8
	sn3d[shi(5, 0)] = z*(63*z**4 - 70*z**2 + 15)/8
	sn3d[shi(5,+1)] = 15**(1./2)*(-z**2 + 1)**(1./2)*(21*z**4 - 14*z**2 + 1)*np.cos(ra)/8
	sn3d[shi(5,+2)] = 105**(1./2)*z*(-3*z**4 + 4*z**2 - 1)*np.cos(2*ra)/4
	sn3d[shi(5,+3)] = 70**(1./2)*(-z**2 + 1)**(1./2)*(-9*z**4 + 10*z**2 - 1)*np.cos(3*ra)/16
	sn3d[shi(5,+4)] = 3*35**(1./2)*z*(z**4 - 2*z**2 + 1)*np.cos(4*ra)/8
	sn3d[shi(5,+5)] = 3*14**(1./2)*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.cos(5*ra)/16
	"""
	sn3d[shi(6,-6)] = 462**(1./2)*(-z**6 + 3*z**4 - 3*z**2 + 1)*np.sin(6*ra)/32
	sn3d[shi(6,-5)] = 3*154**(1./2)*z*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.sin(5*ra)/16
	sn3d[shi(6,-4)] = 3*7**(1./2)*(11*z**6 - 23*z**4 + 13*z**2 - 1)*np.sin(4*ra)/16
	sn3d[shi(6,-3)] = 210**(1./2)*z*(-z**2 + 1)**(1./2)*(-11*z**4 + 14*z**2 - 3)*np.sin(3*ra)/16
	sn3d[shi(6,-2)] = 210**(1./2)*(-33*z**6 + 51*z**4 - 19*z**2 + 1)*np.sin(2*ra)/32
	sn3d[shi(6,-1)] = 21**(1./2)*z*(-z**2 + 1)**(1./2)*(33*z**4 - 30*z**2 + 5)*np.sin(ra)/8
	sn3d[shi(6, 0)] = 231*z**6/16 - 315*z**4/16 + 105*z**2/16 - 5./16
	sn3d[shi(6,+1)] = 21**(1./2)*z*(-z**2 + 1)**(1./2)*(33*z**4 - 30*z**2 + 5)*np.cos(ra)/8
	sn3d[shi(6,+2)] = 210**(1./2)*(-33*z**6 + 51*z**4 - 19*z**2 + 1)*np.cos(2*ra)/32
	sn3d[shi(6,+3)] = 210**(1./2)*z*(-z**2 + 1)**(1./2)*(-11*z**4 + 14*z**2 - 3)*np.cos(3*ra)/16
	sn3d[shi(6,+4)] = 3*7**(1./2)*(11*z**6 - 23*z**4 + 13*z**2 - 1)*np.cos(4*ra)/16
	sn3d[shi(6,+5)] = 3*154**(1./2)*z*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.cos(5*ra)/16
	sn3d[shi(6,+6)] = 462**(1./2)*(-z**6 + 3*z**4 - 3*z**2 + 1)*np.cos(6*ra)/32
	"""

	return sn3d

def semiNormalizedSH_cartesian(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
	x,y,z = ead2xyz(e, a, 1)

	sn3d = np.zeros((order+1,order+1)) if target is None else target
	from sympy import factorial, sqrt, S
	import sympy as sp
	ra = math.radians(a)
	# Synthetic version

	sn3d[shi(0, 0)] = 1.

	sn3d[shi(1,+1)] = x
	sn3d[shi(1, 0)] = z
	sn3d[shi(1,-1)] = y
	sn3d[shi(2,-2)] = (x*y)*2       * np.sqrt(3./4)
	sn3d[shi(2,-1)] = (y*z)*2       * np.sqrt(3./4)
	sn3d[shi(2, 0)] = (3*z*z -1)/2
	sn3d[shi(2,+1)] = (x*z)*2       * np.sqrt(3./4)
	sn3d[shi(2,+2)] = (x*x - y*y)   * np.sqrt(3./4)

	sn3d[shi(3,+3)] = x*(x*x-3*y*y) * np.sqrt(5./8)
	sn3d[shi(3,+2)] = z*(x*x-y*y)   * np.sqrt(15./4)
	sn3d[shi(3,+1)] = x*(5*z*z -1)  * np.sqrt(3./8)
	sn3d[shi(3, 0)] = z*(5*z*z -3)  * np.sqrt(1./4)
	sn3d[shi(3,-1)] = y*(5*z*z -1)  * np.sqrt(3./8)
	sn3d[shi(3,-2)] = z*x*y*2       * np.sqrt(15./4)
	sn3d[shi(3,-3)] = y*(3*x*x-y*y) * np.sqrt(5./8)


#	sn3d[shi(4,+4)] = (x*x*(x*x-3*y*y)-y*y*(3*x*x-y*y))   * np.sqrt(35./64) # reference
	sn3d[shi(4,+4)] = (x**4 -6*x*x*y*y +y**4)  * np.sqrt(35./64)
	sn3d[shi(4,+3)] = z*x*(x*x -3*y*y)         * np.sqrt(35./8)
	sn3d[shi(4,+2)] = (x*x-y*y)*(7*z*z -1)     * np.sqrt(5./16)
	sn3d[shi(4,+1)] = z*x*(7*z*z -3)           * np.sqrt(5./8)
	sn3d[shi(4, 0)] = (3+z*z*(-30+z*z*35))     * np.sqrt(1./64)
	sn3d[shi(4,-1)] = z*y*(7*z*z -3)           * np.sqrt(5./8)
	sn3d[shi(4,-2)] = 2*x*y*(7*z*z -1)         * np.sqrt(5./16)
	sn3d[shi(4,-3)] = z*y*(3*x*x -y*y)         * np.sqrt(35./8)
	sn3d[shi(4,-4)] = 4*x*y*(x*x -y*y)         * np.sqrt(35./64)

	sn3d[shi(5,+5)] = x*(x**4 -10*x*x*y*y +5*y**4)    * np.sqrt(9.*7./64/2)
	sn3d[shi(5,+4)] = z*(x**4 -6*x*x*y*y +y**4)       * np.sqrt(35*9./64)
	sn3d[shi(5,+3)] = x*(x*x-3*y*y)*(9*z*z-1)         * np.sqrt(5.*7./128)
	sn3d[shi(5,+2)] = 3*z*(x*x-y*y)*(3*z*z-1)         * np.sqrt(7.*5/3/16)
	sn3d[shi(5,+1)] = x*(21*z**4 -14*z**2 + 1)        * np.sqrt(3.*5/64)
	sn3d[shi(5, 0)] = z*(63*z**4 -70*z**2 + 15)       * np.sqrt(1./64)
	sn3d[shi(5,-1)] = y*(21*z**4 -14*z**2 + 1)        * np.sqrt(3.*5/64)
	sn3d[shi(5,-2)] = 2*x*y*z*(3*z*z -1)              * np.sqrt(105./16)
	sn3d[shi(5,-3)] = y*(9*z*z-1)*(3*x*x-y*y)         * np.sqrt(35./8/16)
	sn3d[shi(5,-4)] = 4*x*y*z*(x*x -y*y)              * np.sqrt(9.*35./64)
	sn3d[shi(5,-5)] = y*( 5*x**4 -10*x*x*y*y +y**4 )  * np.sqrt(9.*7/64/2)

	assert abs(sn3d[shi(5,+5)]  - 3*14**(1./2)*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.cos(5*ra)/16) < 1e-7
	assert abs(sn3d[shi(5,+4)]  - 3*35**(1./2)*z*(z**4 - 2*z**2 + 1)*np.cos(4*ra)/8 ) < 1e-7
	assert abs(sn3d[shi(5,+3)]  - 70**(1./2)*(-z**2 + 1)**(1./2)*(-9*z**4 + 10*z**2 - 1)*np.cos(3*ra)/16 ) < 1e-7
	assert abs(sn3d[shi(5,+2)]  - 105**(1./2)*z*(-3*z**4 + 4*z**2 - 1)*np.cos(2*ra)/4  ) < 1e-7
	assert abs(sn3d[shi(5,+1)]  - 15**(1./2)*(-z**2 + 1)**(1./2)*(21*z**4 - 14*z**2 + 1)*np.cos(ra)/8  ) < 1e-7
	assert abs(sn3d[shi(5, 0)]  -  z*(63*z**4 - 70*z**2 + 15)/8 )< 1e-7
	assert abs(sn3d[shi(5,-1)]  - 15**(1./2)*(-z**2 + 1)**(1./2)*(21*z**4 - 14*z**2 + 1)*np.sin(ra)/8 ) < 1e-7
	assert abs(sn3d[shi(5,-5)]  - 3*14**(1./2)*(-z**2 + 1)**(1./2)*(z**4 - 2*z**2 + 1)*np.sin(5*ra)/16 ) < 1e-7
	assert abs(sn3d[shi(5,-4)]  - 3*35**(1./2)*z*(z**4 - 2*z**2 + 1)*np.sin(4*ra)/8 ) < 1e-7
	assert abs(sn3d[shi(5,-3)]  - 70**(1./2)*(-z**2 + 1)**(1./2)*(-9*z**4 + 10*z**2 - 1)*np.sin(3*ra)/16 ) < 1e-7
	assert abs(sn3d[shi(5,-2)]  - 105**(1./2)*z*(-3*z**4 + 4*z**2 - 1)*np.sin(2*ra)/4 ) < 1e-7
	assert abs(sn3d[shi(5,-1)]  - 15**(1./2)*(-z**2 + 1)**(1./2)*(21*z**4 - 14*z**2 + 1)*np.sin(ra)/8 ) < 1e-7

	return sn3d

def semiNormalizedSH_cartesianRecursive(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
	x,y,z = ead2xyz(e, a, 1)

	sn3d = np.zeros((order+1,order+1)) if target is None else target
	from sympy import factorial, sqrt, S
	import sympy as sp
	ra = math.radians(a)

	sn3d[shi(0, 0)] = 1.

	sn3d[shi(1,-1)] = y
	sn3d[shi(1, 0)] = z
	sn3d[shi(1,+1)] = x

	sn3d[shi(2,-2)] = (y*sn3d[shi(1,+1)] +x*sn3d[shi(1,-1)]) * np.sqrt(3./4)
	sn3d[shi(2,-1)] = z*sn3d[shi(1,-1)]                      * np.sqrt(3.)
	sn3d[shi(2, 0)] = (3*z*z -1)*sn3d[shi(0, 0)]             * np.sqrt(1./4)
	sn3d[shi(2,+1)] = z*sn3d[shi(1,+1)]                      * np.sqrt(3.)
	sn3d[shi(2,+2)] = (x*sn3d[shi(1,+1)] -y*sn3d[shi(1,-1)]) * np.sqrt(3./4)

	sn3d[shi(3,+3)] = (x*sn3d[shi(2,+2)] -y*sn3d[shi(2,-2)]) * np.sqrt(5./6)
	sn3d[shi(3,+2)] = z*sn3d[shi(2,+2)]                      * np.sqrt(5.)
	sn3d[shi(3,+1)] = (5*z*z -1)*sn3d[shi(1,+1)]             * np.sqrt(3./8)
	sn3d[shi(3, 0)] = (5*z*z -3)*sn3d[shi(1, 0)]             * np.sqrt(1./4)
	sn3d[shi(3,-1)] = (5*z*z -1)*sn3d[shi(1,-1)]             * np.sqrt(3./8)
	sn3d[shi(3,-2)] = z*sn3d[shi(2,-2)]                      * np.sqrt(5.)
	sn3d[shi(3,-3)] = (y*sn3d[shi(2,+2)] +x*sn3d[shi(2,-2)]) * np.sqrt(5./6)

	sn3d[shi(4,+4)] = (x*sn3d[shi(3,+3)] -y*sn3d[shi(3,-3)]) * np.sqrt(7./8)
	sn3d[shi(4,+3)] = z*sn3d[shi(3,+3)]                      * np.sqrt(7.)
	sn3d[shi(4,+2)] = (7*z*z-1)*sn3d[shi(2,+2)]              * np.sqrt(5./3/4)
	sn3d[shi(4,+1)] = (7*z*z-3)*sn3d[shi(2,+1)]              * np.sqrt(5./3/8)
	sn3d[shi(4, 0)] = (3+z*z*(-30+z*z*35))                   * np.sqrt(1./64) # Recursion?
	sn3d[shi(4,-1)] = (7*z*z-3)*sn3d[shi(2,-1)]              * np.sqrt(5./3/8)
	sn3d[shi(4,-2)] = (7*z*z-1)*sn3d[shi(2,-2)]              * np.sqrt(5./3/4)
	sn3d[shi(4,-3)] = z*sn3d[shi(3,-3)]                      * np.sqrt(7.)
	sn3d[shi(4,-4)] = (y*sn3d[shi(3,+3)] +x*sn3d[shi(3,-3)]) * np.sqrt(7./8)

	sn3d[shi(5,+5)] = (x*sn3d[shi(4,+4)] -y*sn3d[shi(4,-4)]) * np.sqrt(9./10)
	sn3d[shi(5,+4)] = z*sn3d[shi(4,+4)]                      * np.sqrt(9.)
	sn3d[shi(5,+3)] = (9*z*z-1)*sn3d[shi(3,+3)]              * np.sqrt(7./8/2)
	sn3d[shi(5,+2)] = (9*z*z-3)*sn3d[shi(3,+2)]              * np.sqrt(7./9/4)
	sn3d[shi(5,+1)] = x*(21*z**4 -14*z**2 + 1)               * np.sqrt(3.*5/64)
	sn3d[shi(5, 0)] = z*(63*z**4 -70*z**2 + 15)              * np.sqrt(1./64)
	sn3d[shi(5,-1)] = y*(21*z**4 -14*z**2 + 1)               * np.sqrt(3.*5/64)
	sn3d[shi(5,-2)] = 2*x*y*z*(3*z*z -1)                     * np.sqrt(105./16)
	sn3d[shi(5,-3)] = y*(9*z*z-1)*(3*x*x-y*y)                * np.sqrt(35./8/16)
	sn3d[shi(5,-4)] = z*sn3d[shi(4,-4)]                      * np.sqrt(9.)
	sn3d[shi(5,-5)] = (x*sn3d[shi(4,-4)] +y*sn3d[shi(4,+4)]) * np.sqrt(9./10)

	return sn3d


def semiNormalizedSH(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
	return semiNormalizedSH_cartesian(e,a,target)

def sh(components, e, a) :
	"""
	Decodes the sh components for (e,a) position.
	"""
	projection = semiNormalizedSH(e,a)
	return (components*projection).sum()



if __name__ == "__main__" :
	unittest.main()
	width = 600
	height = 400

	app = QtGui.QApplication(sys.argv)

	w = QtGui.QDialog()
	w.setLayout(QtGui.QVBoxLayout())
	w1 = ColorField(width, height)
	reloader1 = Reloader(w1)
	reloader1.startTimer(0)
	w.layout().addWidget(w1)

	w2 = ColorField(width, height, ColorField.fancyScale)
	reloader2 = Reloader(w2)
	reloader2.startTimer(0)
	w.layout().addWidget(w2)

	w.show()
	w.resize(width, 2*height)

	app.exec_()

