#!/usr/bin/python

from PySide import QtCore, QtGui, QtOpenGL
import sys
import numpy as np
import time
import math

order = 5
shShape = (order+1,order+1)
shSize = np.prod(shShape)

def shIndexes(order) :
	return [(l,m) for l in xrange(order+1) for m in xrange(-l,l+1) ]

def shi(l,m) :
	return (l-m,l) if m>0 else (l,l+m)

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


class SemiNormalizedSH(object) :
	shXYZ = dict()
	def __init__(self, compiled=False) :
		from sympy import S, symbols, sqrt, Lambda
		from sympy.utilities.autowrap import autowrap
		x,y,z = symbols("x y z")

		self.order = order
		self.sn3d = np.array([[None]*(self.order+1)]*(self.order+1))

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


		if compiled :
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
					self.sn3d[i,j] = Lambda(
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

_associated_legendre_cache = {}
def associated_legendre(l,m,z) :
	import sympy as sp
	if (l,m,z) in _associated_legendre_cache :
		print "using cache"
		return _associated_legendre_cache[(l,m,z)]
	return sp.assoc_legendre(l,m,z)

def semiNormalizedSH(e, a, target=None) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
#	return obj.evalEA(e,a, target)

	x,y,z = ead2xyz(e, a, 1)

	sn3d = np.zeros((order+1,order+1)) if target is None else target
	from sympy import factorial, sqrt, S
	import sympy as sp
	ra = math.radians(a)
	"""
	for l in xrange(order+1) :
		for m in xrange(-l,l+1) :
			absm = abs(m)
			factor = math.sqrt(2) if m else 1.
			sn3d[shi(l,m)] = (factor
				/ math.sqrt(
					sp.factorial(l+abs(m)) / sp.factorial(l-abs(m)) # clearest for the next one
#					np.prod(xrange(l-absm+1, l+absm+1))
					)
				* ( sp.assoc_legendre(l,absm,z) * S(-1)**m ).evalf()
				* ( np.cos(m*ra) if m>=0 else np.sin(-m*ra) )
				)
#	return sn3d
	"""

	"""
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

	"""
# recursive versions

#	sn3d[shi(1,-1)] = (y*sn3d[shi(0,+0)] +x*sn3d[shi(0,-0)]) * np.sqrt(1)
#	sn3d[shi(1, 0)] = (z*z -1)
#	sn3d[shi(1,+1)] = (x*sn3d[shi(0,+0)] -y*sn3d[shi(0,-0)]) * np.sqrt(1)

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

	"""

	return sn3d[:order+1,:order+1]


def sh(components, e, a) :
	"""
	Decodes the sh components for (e,a) position.
	"""
	projection = semiNormalizedSH(e,a)
	return (components*projection).sum()

########################################
# Normalization factors,
# They are to be applied to the sn3d to get a different one.


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

# TODO: maxn for 4th order
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



import unittest


class SphericalHarmonicsTests(unittest.TestCase) :

	def _fumaNormalizedComponents(self, l, m) :
		"""
		Returns a SH projection that has only the indicated component
		and is normalized using the FuMa normalization so that the maximum
		value is one (but the 0,0 component which is normalized to 1/sqrt(2))
		"""
		components = np.zeros(shShape)
		components[shi(l,m)] = fuma[shi(l,m)]
		self._components = components
		return components

	def prepareOrder(self, l,m) :
		self._components = self._fumaNormalizedComponents(l,m)

	def assertEqualAt(self, ea, value) :
		e,a=ea
		self.assertAlmostEqual(sh(self._components, e, a), value)

	def test_shIndexes(self) :
		self.assertEqual(
			[
				(0,0),
				(1,-1), (1,0), (1,1),
				(2,-2), (2,-1), (2,0), (2,1), (2,2)
			],
			shIndexes(2),
			)

	def test_shIndex2Matrix(self) :

		self.assertEqual((0,0), shi(0, 0))

		self.assertEqual((1,0), shi(1,-1))
		self.assertEqual((1,1), shi(1, 0))
		self.assertEqual((0,1), shi(1,+1))

		self.assertEqual((2,0), shi(2,-2))
		self.assertEqual((2,1), shi(2,-1))
		self.assertEqual((2,2), shi(2, 0))
		self.assertEqual((1,2), shi(2,+1))
		self.assertEqual((0,2), shi(2,+2))


	def test_sh_0_0(self) :
		self.prepareOrder(0, 0)

		maxValue = math.sqrt(1./2) # FuMa is not MaxN for 0,0

		self.assertEqualAt((+90, 0), +maxValue)
		self.assertEqualAt((-90, 0), +maxValue)

		self.assertEqualAt(( 30, 4), +maxValue)

	def test_sh_1_1(self) :
		self.prepareOrder(1,+1)

		self.assertEqualAt((   0,   0),+1)
		self.assertEqualAt((   0, 180),-1)
		self.assertEqualAt((   0, -90), 0)
		self.assertEqualAt((   0, +90), 0)
		self.assertEqualAt(( -90,   0), 0)
		self.assertEqualAt(( +90,   0), 0)
		self.assertEqualAt((  30,   4), 0.86391580942710433)

	def test_sh_1_m1(self) :
		self.prepareOrder(1,-1)

		self.assertEqualAt((   0,   0), 0)
		self.assertEqualAt((   0, 180), 0)
		self.assertEqualAt((   0, -90),-1)
		self.assertEqualAt((   0, +90),+1)
		self.assertEqualAt(( -90,   0), 0)
		self.assertEqualAt(( +90,   0), 0)
		self.assertEqualAt((  30,   4), 0.060410878340834702)

	def test_sh_1_0(self) :
		self.prepareOrder(1,0)

		self.assertEqualAt((   0,   0), 0)
		self.assertEqualAt((   0, 180), 0)
		self.assertEqualAt((   0, -90), 0)
		self.assertEqualAt((   0, +90), 0)
		self.assertEqualAt(( -90,   0),-1)
		self.assertEqualAt(( +90,   0),+1)
		self.assertEqualAt((  30,   0), 0.5)

		self.assertEqualAt((  30,   4), 0.49999999999999994)

	def test_sh_2_2(self) :
		self.prepareOrder(2,+2)

		self.assertEqualAt((   0,   0), +1)
		self.assertEqualAt((   0, 180), +1)
		self.assertEqualAt((   0, -90), -1)
		self.assertEqualAt((   0, +90), -1)
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		self.assertEqualAt((  45,   0), .5)
		self.assertEqualAt((  45, 180), .5)
		self.assertEqualAt(( -45,   0), .5)
		self.assertEqualAt(( -45, 180), .5)

		self.assertEqualAt((   0,  45),  0)
		self.assertEqualAt((   0, 135),  0)
		self.assertEqualAt((   0, 225),  0)
		self.assertEqualAt((   0, 315),  0)

		self.assertEqualAt((  45, -90), -.5)
		self.assertEqualAt((  45, +90), -.5)
		self.assertEqualAt(( -45, -90), -.5)
		self.assertEqualAt(( -45, +90), -.5)

		self.assertEqualAt((  30,   4), 0.7427010515561776)

	def test_sh_2_1(self) :
		self.prepareOrder(2,+1)

		# axis
		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)
		self.assertEqualAt((   0, -90),  0)
		self.assertEqualAt((   0, +90),  0)
		# extremes
		self.assertEqualAt((  45,   0), +1)
		self.assertEqualAt((  45, 180), -1)
		self.assertEqualAt(( -45,   0), -1)
		self.assertEqualAt(( -45, 180), +1)

		self.assertEqualAt((   0,  45),  0)
		self.assertEqualAt((   0, 135),  0)
		self.assertEqualAt((   0, 225),  0)
		self.assertEqualAt((   0, 315),  0)

		self.assertEqualAt((  45, -90),  0)
		self.assertEqualAt((  45, +90),  0)
		self.assertEqualAt(( -45, -90),  0)
		self.assertEqualAt(( -45, +90),  0)

		self.assertEqualAt((  30,   4),  0.86391580942710411)

	def test_sh_2_0(self) :
		self.prepareOrder(2, 0)

		# axis
		self.assertEqualAt((   0,   0), -.5)
		self.assertEqualAt((   0, 180), -.5)
		self.assertEqualAt((   0, -90), -.5)
		self.assertEqualAt((   0, +90), -.5)
		self.assertEqualAt(( -90,   0),  +1)
		self.assertEqualAt(( +90,   0),  +1)
		# extremes
		self.assertEqualAt((   0,  45), -.5)
		self.assertEqualAt((   0, 135), -.5)
		self.assertEqualAt((   0, 225), -.5)
		self.assertEqualAt((   0, 315), -.5)

		self.assertEqualAt((  45, -90), +.25)
		self.assertEqualAt((  45, +90), +.25)
		self.assertEqualAt(( -45, -90), +.25)
		self.assertEqualAt(( -45, +90), +.25)

		self.assertEqualAt((  45,   0), +.25)
		self.assertEqualAt((  45, 180), +.25)
		self.assertEqualAt(( -45,   0), +.25)
		self.assertEqualAt(( -45, 180), +.25)

		self.assertEqualAt((  30,   4),  -0.1250000000000011)

	def test_sh_2_m1(self) :
		self.prepareOrder(2,-1)

		# axis
		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)
		self.assertEqualAt((   0, -90),  0)
		self.assertEqualAt((   0, +90),  0)
		# extremes
		self.assertEqualAt((   0,  45),  0)
		self.assertEqualAt((   0, 135),  0)
		self.assertEqualAt((   0, 225),  0)
		self.assertEqualAt((   0, 315),  0)

		self.assertEqualAt((  45, -90), -1)
		self.assertEqualAt((  45, +90), +1)
		self.assertEqualAt(( -45, -90), +1)
		self.assertEqualAt(( -45, +90), -1)

		self.assertEqualAt((  45,   0),  0)
		self.assertEqualAt((  45, 180),  0)
		self.assertEqualAt(( -45,   0),  0)
		self.assertEqualAt(( -45, 180),  0)

		self.assertEqualAt((  30,   4),  0.060410878340834695)


	def test_sh_2_m2(self) :
		self.prepareOrder(2,-2)

		# axis
		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)
		self.assertEqualAt((   0, -90),  0)
		self.assertEqualAt((   0, +90),  0)
		# extremes
		self.assertEqualAt((   0,  45), +1)
		self.assertEqualAt((   0, 135), -1)
		self.assertEqualAt((   0, 225), +1)
		self.assertEqualAt((   0, 315), -1)

		self.assertEqualAt((  45, -90),  0)
		self.assertEqualAt((  45, +90),  0)
		self.assertEqualAt(( -45, -90),  0)
		self.assertEqualAt(( -45, +90),  0)

		self.assertEqualAt((  45,   0),  0)
		self.assertEqualAt((  45, 180),  0)
		self.assertEqualAt(( -45,   0),  0)
		self.assertEqualAt(( -45, 180),  0)

		self.assertEqualAt((  30,   4),  0.10437982572004907)

	def test_sh_3_3(self) :
		self.prepareOrder(3,+3)

		# axis
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)
		self.assertEqualAt((   0,   0), +1)
		self.assertEqualAt((   0,  30),  0)
		self.assertEqualAt((   0,  60), -1)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 120), +1)
		self.assertEqualAt((   0, 150),  0)
		self.assertEqualAt((   0, 180), -1)
		self.assertEqualAt((   0, 210),  0)
		self.assertEqualAt((   0, 240), +1)
		self.assertEqualAt((   0, 270),  0)
		self.assertEqualAt((   0, 300), -1)
		self.assertEqualAt((   0, 330),  0)

		self.assertEqualAt((  30,   4),  0.63532550316470549)

	def test_sh_3_2(self) :
		self.prepareOrder(3,+2)

		# axis
		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)
		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0,  30),  0)
		self.assertEqualAt((   0,  60),  0)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 120),  0)
		self.assertEqualAt((   0, 150),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt((   0, 210),  0)
		self.assertEqualAt((   0, 240),  0)
		self.assertEqualAt((   0, 270),  0)
		self.assertEqualAt((   0, 300),  0)
		self.assertEqualAt((   0, 330),  0)

		maxangle = math.degrees(math.asin(math.sqrt(1./3))) # 35.2644

		self.assertEqualAt(( +maxangle,   0), +1)
		self.assertEqualAt(( +maxangle,  45),  0)
		self.assertEqualAt(( +maxangle,  90), -1)
		self.assertEqualAt(( +maxangle, 135),  0)
		self.assertEqualAt(( +maxangle, 180), +1)
		self.assertEqualAt(( +maxangle, 225),  0)
		self.assertEqualAt(( +maxangle, 270), -1)
		self.assertEqualAt(( +maxangle, 315),  0)

		self.assertEqualAt(( -maxangle,   0), -1)
		self.assertEqualAt(( -maxangle,  45),  0)
		self.assertEqualAt(( -maxangle,  90), +1)
		self.assertEqualAt(( -maxangle, 135),  0)
		self.assertEqualAt(( -maxangle, 180), -1)
		self.assertEqualAt(( -maxangle, 225),  0)
		self.assertEqualAt(( -maxangle, 270), +1)
		self.assertEqualAt(( -maxangle, 315),  0)

		self.assertEqualAt((  30,   4),  0.96479696709759877)

	def test_sh_3_1(self) :
		self.prepareOrder(3,+1)

		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		xymaxvalue = math.sqrt(3*5)*3/16 # -0.72618437741389052

		self.assertEqualAt((   0,   0), -xymaxvalue)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 180), +xymaxvalue)
		self.assertEqualAt((   0, 270),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(1./5)))

		self.assertEqualAt(( +zeroangle,   0),  0)
		self.assertEqualAt(( +zeroangle,  90),  0)
		self.assertEqualAt(( +zeroangle, 180),  0)
		self.assertEqualAt(( +zeroangle, 270),  0)

		self.assertEqualAt(( -zeroangle,   0),  0)
		self.assertEqualAt(( -zeroangle,  90),  0)
		self.assertEqualAt(( -zeroangle, 180),  0)
		self.assertEqualAt(( -zeroangle, 270),  0)

		maxangle = math.degrees(math.acos(math.sqrt(4./15))) # 58.90907

		self.assertEqualAt(( +maxangle,   0), +1)
		self.assertEqualAt(( +maxangle,  90),  0)
		self.assertEqualAt(( +maxangle, 180), -1)
		self.assertEqualAt(( +maxangle, 270),  0)

		self.assertEqualAt(( -maxangle,   0), +1)
		self.assertEqualAt(( -maxangle,  90),  0)
		self.assertEqualAt(( -maxangle, 180), -1)
		self.assertEqualAt(( -maxangle, 270),  0)

		self.assertEqualAt((  30,   4),  0.15684054105170947)

	def test_sh_3_0(self) :
		self.prepareOrder(3, 0)

		# axis
		self.assertEqualAt(( -90,   0), -1)
		self.assertEqualAt(( +90,   0), +1)
		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0,  30),  0)
		self.assertEqualAt((   0,  60),  0)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 120),  0)
		self.assertEqualAt((   0, 150),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt((   0, 210),  0)
		self.assertEqualAt((   0, 240),  0)
		self.assertEqualAt((   0, 270),  0)
		self.assertEqualAt((   0, 300),  0)
		self.assertEqualAt((   0, 330),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(3./5))) # 50.7684795164

		self.assertEqualAt(( +zeroangle,   0),  0)
		self.assertEqualAt(( +zeroangle,  90),  0)
		self.assertEqualAt(( +zeroangle, 180),  0)
		self.assertEqualAt(( +zeroangle, 270),  0)

		self.assertEqualAt(( -zeroangle,   0),  0)
		self.assertEqualAt(( -zeroangle,  90),  0)
		self.assertEqualAt(( -zeroangle, 180),  0)
		self.assertEqualAt(( -zeroangle, 270),  0)

		maxangle = math.degrees(math.asin(math.sqrt(3./15))) # 26.5650511771
		maxvalue = math.sqrt(3./15)

		self.assertEqualAt(( +maxangle,   0), -maxvalue)
		self.assertEqualAt(( +maxangle,  90), -maxvalue)
		self.assertEqualAt(( +maxangle, 180), -maxvalue)
		self.assertEqualAt(( +maxangle, 270), -maxvalue)

		self.assertEqualAt(( -maxangle,   0), +maxvalue)
		self.assertEqualAt(( -maxangle,  90), +maxvalue)
		self.assertEqualAt(( -maxangle, 180), +maxvalue)
		self.assertEqualAt(( -maxangle, 270), +maxvalue)

		self.assertEqualAt((  30,   4),  -0.43750000000000006)

	def test_sh_3_m1(self) :
		self.prepareOrder(3,-1)

		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		xymaxvalue = math.sqrt(3*5)*3/16 # -0.72618437741389052
		self.assertEqualAt((   0,  90), -xymaxvalue)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt((   0, 270), +xymaxvalue)
		self.assertEqualAt((   0,   0),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(1./5)))

		self.assertEqualAt(( +zeroangle,   0),  0)
		self.assertEqualAt(( +zeroangle,  90),  0)
		self.assertEqualAt(( +zeroangle, 180),  0)
		self.assertEqualAt(( +zeroangle, 270),  0)

		self.assertEqualAt(( -zeroangle,   0),  0)
		self.assertEqualAt(( -zeroangle,  90),  0)
		self.assertEqualAt(( -zeroangle, 180),  0)
		self.assertEqualAt(( -zeroangle, 270),  0)

		maxangle = math.degrees(math.acos(math.sqrt(4./15))) # 58.90907
		self.assertEqualAt(( +maxangle,   0),  0)
		self.assertEqualAt(( +maxangle,  90), +1)
		self.assertEqualAt(( +maxangle, 180),  0)
		self.assertEqualAt(( +maxangle, 270), -1)

		self.assertEqualAt(( -maxangle,   0),  0)
		self.assertEqualAt(( -maxangle,  90), +1)
		self.assertEqualAt(( -maxangle, 180),  0)
		self.assertEqualAt(( -maxangle, 270), -1)

		self.assertEqualAt((  30,   4),  0.010967359019241315)

	def test_sh_3_m2(self) :
		self.prepareOrder(3,-2)

		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0,  30),  0)
		self.assertEqualAt((   0,  60),  0)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 120),  0)
		self.assertEqualAt((   0, 150),  0)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt((   0, 210),  0)
		self.assertEqualAt((   0, 240),  0)
		self.assertEqualAt((   0, 270),  0)
		self.assertEqualAt((   0, 300),  0)
		self.assertEqualAt((   0, 330),  0)


		maxangle = math.degrees(math.asin(math.sqrt(1./3))) # 35.2644

		self.assertEqualAt(( +maxangle,   0),  0)
		self.assertEqualAt(( +maxangle,  45), +1)
		self.assertEqualAt(( +maxangle,  90),  0)
		self.assertEqualAt(( +maxangle, 135), -1)
		self.assertEqualAt(( +maxangle, 180),  0)
		self.assertEqualAt(( +maxangle, 225), +1)
		self.assertEqualAt(( +maxangle, 270),  0)
		self.assertEqualAt(( +maxangle, 315), -1)

		self.assertEqualAt(( -maxangle,   0),  0)
		self.assertEqualAt(( -maxangle,  45), -1)
		self.assertEqualAt(( -maxangle,  90),  0)
		self.assertEqualAt(( -maxangle, 135), +1)
		self.assertEqualAt(( -maxangle, 180),  0)
		self.assertEqualAt(( -maxangle, 225), -1)
		self.assertEqualAt(( -maxangle, 270),  0)
		self.assertEqualAt(( -maxangle, 315), +1)

		self.assertEqualAt((  30,   4),  0.13559337107423225)

	def test_sh_3_m3(self) :
		self.prepareOrder(3,-3)

		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		self.assertEqualAt((   0,   0),  0)
		self.assertEqualAt((   0,  30), +1)
		self.assertEqualAt((   0,  60),  0)
		self.assertEqualAt((   0,  90), -1)
		self.assertEqualAt((   0, 120),  0)
		self.assertEqualAt((   0, 150), +1)
		self.assertEqualAt((   0, 180),  0)
		self.assertEqualAt((   0, 210), -1)
		self.assertEqualAt((   0, 240),  0)
		self.assertEqualAt((   0, 270), +1)
		self.assertEqualAt((   0, 300),  0)
		self.assertEqualAt((   0, 330), -1)

		self.assertEqualAt((  30,   4),  0.13504260449396654)

	# TODO: Warning 4th and above are not FuMa normalized yet!!!

	def test_sh_4_0(self) :
		self.prepareOrder(4,-0)

		self.assertEqualAt(( -90,   0), +1)
		self.assertEqualAt(( +90,   0), +1)

		maxvalue = 3./8

		self.assertEqualAt(( 0,   0), maxvalue)
		self.assertEqualAt(( 0,  90), maxvalue)
		self.assertEqualAt(( 0, 180), maxvalue)
		self.assertEqualAt(( 0, -90), maxvalue)

		maxangle = math.degrees(math.asin(math.sqrt(3./7))) # 40.8933946491309
		maxvalue = -3./7

		self.assertEqualAt(( maxangle,   0), maxvalue)
		self.assertEqualAt(( maxangle,  90), maxvalue)
		self.assertEqualAt(( maxangle, 180), maxvalue)
		self.assertEqualAt(( maxangle, -90), maxvalue)

		zero1 = math.degrees(math.asin(math.sqrt(3./7 + 2*math.sqrt(6./5)/7)))
		zero2 = math.degrees(math.asin(math.sqrt(3./7 - 2*math.sqrt(6./5)/7)))

		self.assertEqualAt(( -zero1,   0),  0)
		self.assertEqualAt(( -zero2,   0),  0)
		self.assertEqualAt(( +zero2,   0),  0)
		self.assertEqualAt(( +zero1,   0),  0)

		self.assertEqualAt((  30,   4), -0.28906249999999989)

	def test_sh_4_1(self) :
		self.prepareOrder(4,+1)

		self.assertEqualAt((  30,   4), -0.42686588506525242)

		self.assertEqualAt(( -90,   0), 0)
		self.assertEqualAt(( +90,   0), 0)

		self.assertEqualAt(( 0,   0), 0)
		self.assertEqualAt(( 0,  30), 0)
		self.assertEqualAt(( 0,  45), 0)
		self.assertEqualAt(( 0,  60), 0)
		self.assertEqualAt(( 0,  90), 0)
		self.assertEqualAt(( 0, 180), 0)
		self.assertEqualAt(( 0, -90), 0)

		zeroangle = math.degrees(math.asin(math.sqrt(3./7)))

		self.assertEqualAt(( -zeroangle, +90),  0)
		self.assertEqualAt(( -zeroangle, + 3),  0)
		self.assertEqualAt(( +zeroangle, +90),  0)
		self.assertEqualAt(( +zeroangle, + 3),  0)

		maxangle1 = math.degrees(math.asin(math.sqrt( (27-math.sqrt(393))/7/8) ))
		maxvalue1 = 0.55571119769376354 # b2b, not analytically checked
		maxangle2 = math.degrees(math.asin(math.sqrt( (27+math.sqrt(393))/7/8) )) # 66.1221762567095 
		maxvalue2 = 0.83486203359622035 # b2b, not analytically checked

		self.assertEqualAt(( +maxangle1,   0), -maxvalue1)
		self.assertEqualAt(( +maxangle1,  90),  0)
		self.assertEqualAt(( +maxangle1, 180), +maxvalue1)
		self.assertEqualAt(( +maxangle1, -90),  0)

		self.assertEqualAt(( +maxangle2,   0), +maxvalue2)
		self.assertEqualAt(( +maxangle2,  90),  0)
		self.assertEqualAt(( +maxangle2, 180), -maxvalue2)
		self.assertEqualAt(( +maxangle2, -90),  0)

		self.assertEqualAt(( -maxangle1,   0), +maxvalue1)
		self.assertEqualAt(( -maxangle1,  90),  0)
		self.assertEqualAt(( -maxangle1, 180), -maxvalue1)
		self.assertEqualAt(( -maxangle1, -90),  0)

		self.assertEqualAt(( -maxangle2,   0), -maxvalue2)
		self.assertEqualAt(( -maxangle2,  90),  0)
		self.assertEqualAt(( -maxangle2, 180), +maxvalue2)
		self.assertEqualAt(( -maxangle2, -90),  0)

	def test_sh_4_2(self) :
		self.prepareOrder(4,+2)

		self.assertEqualAt((  30,   4), 0.3113868821700353)

	def test_sh_4_3(self) :
		self.prepareOrder(4,+3)

		self.assertEqualAt((  30,   4), 0.66443931541944656)

	def test_sh_4_4(self) :
		self.prepareOrder(4,+4)

		self.assertEqualAt((  30,   4), 0.39986021851936454)

	def test_sh_4_m1(self) :
		self.prepareOrder(4,-1)

		self.assertEqualAt((  30,   4), -0.029849370470058041)

	def test_sh_4_m2(self) :
		self.prepareOrder(4,-2)

		self.assertEqualAt((  30,   4), 0.043762572335551975)

	def test_sh_4_m3(self) :
		self.prepareOrder(4,-3)

		self.assertEqualAt((  30,   4), 0.14123093632394088)

	def test_sh_4_m4(self) :
		self.prepareOrder(4,-4)

		self.assertEqualAt((  30,   4), 0.1146580726089364)

	def test_sh_5_all(self) :
		self.prepareOrder(5,-5)
		self.assertEqualAt((  30,   4), +0.116888055250392)
		self.prepareOrder(5,-4)
		self.assertEqualAt((  30,   4), +0.171987108913405)
		self.prepareOrder(5,-3)
		self.assertEqualAt((  30,   4), +0.0882693352024631)
		self.prepareOrder(5,-2)
		self.assertEqualAt((  30,   4), -0.0334242167222746)
		self.prepareOrder(5,-1)
		self.assertEqualAt((  30,   4), -0.0347299702275976)
		self.prepareOrder(5, 0)
		self.assertEqualAt((  30,   4), +0.0898437500000000)
		self.prepareOrder(5,+1)
		self.assertEqualAt((  30,   4), -0.496661713330414)
		self.prepareOrder(5,+2)
		self.assertEqualAt((  30,   4), -0.237825659660080)
		self.prepareOrder(5,+3)
		self.assertEqualAt((  30,   4), +0.415274572137154)
		self.prepareOrder(5,+4)
		self.assertEqualAt((  30,   4), +0.599790327779047)
		self.prepareOrder(5,+5)
		self.assertEqualAt((  30,   4), +0.321147292404416)


class CoordsConversionTests(unittest.TestCase) :

	def assertCoordsEqual(self, expected, result) :
		self.assertEqual(
			[round(v,7) for v in expected],
			[round(v,7) for v in result],
			)

	def test_ead2xyz_front(self) :
		self.assertCoordsEqual( [1.0, 0.0, 0.0], ead2xyz(0, 0, 1) )

	def test_ead2xyz_further_front(self) :
		self.assertCoordsEqual( [2.0, 0.0, 0.0], ead2xyz(0, 0, 2) )

	def test_ead2xyz_back_turning_left(self) :
		self.assertCoordsEqual( [-1.0, 0.0, 0.0], ead2xyz(0, 180, 1) )

	def test_ead2xyz_back_turning_right(self) :
		self.assertCoordsEqual( [-1.0, 0.0, 0.0], ead2xyz(0, -180, 1) )

	def test_ead2xyz_left(self) :
		self.assertCoordsEqual( [0.0, 1.0, 0.0], ead2xyz(0, 90, 1) )

	def test_ead2xyz_right(self) :
		self.assertCoordsEqual( [0.0, -1.0, 0.0], ead2xyz(0, -90, 1) )

	def test_ead2xyz_up(self) :
		self.assertCoordsEqual( [0.0, 0.0, 1.0], ead2xyz(90, 0, 1) )

	def test_ead2xyz_up_despite_azimuth(self) :
		self.assertCoordsEqual( [0.0, 0.0, 1.0], ead2xyz(90, 80, 1) )

	def test_ead2xyz_down(self) :
		self.assertCoordsEqual( [0.0, 0.0, -1.0], ead2xyz(-90, 0, 1) )

	def test_ead2xyz_mid_elevated(self) :
		self.assertCoordsEqual( [0.8660254, 0.0, 0.5], ead2xyz(30, 0, 1) )

	def test_ead2xyz_mixed(self) :
		self.assertCoordsEqual( [0.75, 0.4330127, 0.5], ead2xyz(30, 30, 1) )


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

