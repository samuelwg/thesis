#!/usr/bin/python

from PySide import QtCore, QtGui, QtOpenGL
import sys
import numpy as np
import time
import math

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

	ra, re = math.radians(a), math.radians(e)
	sa, se = math.sin(ra), math.sin(re)
	ca, ce = math.cos(ra), math.cos(re)
	x,y,z = d*ce*ca, d*ce*sa, d*se
	return x,y,z

def semiNormalizedSH(e, a) :
	"""
	Returns the value of the sh components at the specified orientation
	"""
	x,y,z = ead2xyz(e, a, 1)
	sn3d = np.zeros((4,4))

	sn3d[shi(0, 0)] = 1.

	sn3d[shi(1,+1)] = x
	sn3d[shi(1, 0)] = z
	sn3d[shi(1,-1)] = y

	sn3d[shi(2,-2)] = (x*y)*2       * math.sqrt(3./4)
	sn3d[shi(2,-1)] = (y*z)*2       * math.sqrt(3./4)
	sn3d[shi(2, 0)] = (3*z*z -1)/2
	sn3d[shi(2,+1)] = (x*z)*2       * math.sqrt(3./4)
	sn3d[shi(2,+2)] = (x*x - y*y)   * math.sqrt(3./4)

	sn3d[shi(3,+3)] = x*(x*x-3*y*y) * math.sqrt(5./8)
	sn3d[shi(3,+2)] = z*(x*x-y*y)   * math.sqrt(15./4)
	sn3d[shi(3,+1)] = x*(5*z*z -1)  * math.sqrt(3./8)
	sn3d[shi(3, 0)] = z*(5*z*z -3)  * math.sqrt(1./4)
	sn3d[shi(3,-1)] = y*(5*z*z -1)  * math.sqrt(3./8)
	sn3d[shi(3,-2)] = z*x*y         * math.sqrt(15.)
	sn3d[shi(3,-3)] = y*(3*x*x-y*y) * math.sqrt(5./8)

	return sh


def sh(sh, e, a) :
	"""
	Decodes the sh components for (e,a) position.
	"""
	projection = semiNormalizedSH(e,a)
	return (sh*projection).sum()


# factors to be applied to a sn3d sh to get n3d normalization.
# n3d normalization is 
n3d = np.zeros((4,4))
for l in xrange(0,4) :
	for m in xrange(-l, l+1) :
		n3d[shi(l,m)] = math.sqrt(2*l+1)

# factors to be applied to sn3d to get maxn normalization.
# maxn normalization is that so that absolute maximum value
# for each sh is 1.
maxn=np.zeros((4,4))
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
		components = np.zeros((4,4))
		components[shi(l,m)] = fuma[shi(l,m)]
		return components

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
		components = self._fumaNormalizedComponents(0, 0)

		maxValue = math.sqrt(1./2) # FuMa is not MaxN for 0,0

		self.assertAlmostEqual(sh(components,+90, 0), +maxValue)
		self.assertAlmostEqual(sh(components,-90, 0), +maxValue)

		self.assertAlmostEqual(sh(components, 30, 4), +maxValue)

	def test_sh_1_1(self) :
		components = self._fumaNormalizedComponents(1,+1)

		self.assertAlmostEqual(sh(components,   0,   0),+1)
		self.assertAlmostEqual(sh(components,   0, 180),-1)
		self.assertAlmostEqual(sh(components,   0, -90), 0)
		self.assertAlmostEqual(sh(components,   0, +90), 0)
		self.assertAlmostEqual(sh(components, -90,   0), 0)
		self.assertAlmostEqual(sh(components, +90,   0), 0)
		self.assertAlmostEqual(sh(components,  30,   4), 0.86391580942710433)

	def test_sh_1_m1(self) :
		components = self._fumaNormalizedComponents(1,-1)

		self.assertAlmostEqual(sh(components,   0,   0), 0)
		self.assertAlmostEqual(sh(components,   0, 180), 0)
		self.assertAlmostEqual(sh(components,   0, -90),-1)
		self.assertAlmostEqual(sh(components,   0, +90),+1)
		self.assertAlmostEqual(sh(components, -90,   0), 0)
		self.assertAlmostEqual(sh(components, +90,   0), 0)
		self.assertAlmostEqual(sh(components,  30,   4), 0.060410878340834702)

	def test_sh_1_0(self) :
		components = self._fumaNormalizedComponents(1,0)

		self.assertAlmostEqual(sh(components,   0,   0), 0)
		self.assertAlmostEqual(sh(components,   0, 180), 0)
		self.assertAlmostEqual(sh(components,   0, -90), 0)
		self.assertAlmostEqual(sh(components,   0, +90), 0)
		self.assertAlmostEqual(sh(components, -90,   0),-1)
		self.assertAlmostEqual(sh(components, +90,   0),+1)
		self.assertAlmostEqual(sh(components,  30,   0), 0.5)

		self.assertAlmostEqual(sh(components,  30,   4), 0.49999999999999994)

	def test_sh_2_2(self) :
		components = self._fumaNormalizedComponents(2,+2)

		self.assertAlmostEqual(sh(components,   0,   0), +1)
		self.assertAlmostEqual(sh(components,   0, 180), +1)
		self.assertAlmostEqual(sh(components,   0, -90), -1)
		self.assertAlmostEqual(sh(components,   0, +90), -1)
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)

		self.assertAlmostEqual(sh(components,  45,   0), .5)
		self.assertAlmostEqual(sh(components,  45, 180), .5)
		self.assertAlmostEqual(sh(components, -45,   0), .5)
		self.assertAlmostEqual(sh(components, -45, 180), .5)

		self.assertAlmostEqual(sh(components,   0,  45),  0)
		self.assertAlmostEqual(sh(components,   0, 135),  0)
		self.assertAlmostEqual(sh(components,   0, 225),  0)
		self.assertAlmostEqual(sh(components,   0, 315),  0)

		self.assertAlmostEqual(sh(components,  45, -90), -.5)
		self.assertAlmostEqual(sh(components,  45, +90), -.5)
		self.assertAlmostEqual(sh(components, -45, -90), -.5)
		self.assertAlmostEqual(sh(components, -45, +90), -.5)

		self.assertAlmostEqual(sh(components,  30,   4), 0.7427010515561776)

	def test_sh_2_1(self) :
		components = self._fumaNormalizedComponents(2,+1)

		# axis
		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)
		self.assertAlmostEqual(sh(components,   0, -90),  0)
		self.assertAlmostEqual(sh(components,   0, +90),  0)
		# extremes
		self.assertAlmostEqual(sh(components,  45,   0), +1)
		self.assertAlmostEqual(sh(components,  45, 180), -1)
		self.assertAlmostEqual(sh(components, -45,   0), -1)
		self.assertAlmostEqual(sh(components, -45, 180), +1)

		self.assertAlmostEqual(sh(components,   0,  45),  0)
		self.assertAlmostEqual(sh(components,   0, 135),  0)
		self.assertAlmostEqual(sh(components,   0, 225),  0)
		self.assertAlmostEqual(sh(components,   0, 315),  0)

		self.assertAlmostEqual(sh(components,  45, -90),  0)
		self.assertAlmostEqual(sh(components,  45, +90),  0)
		self.assertAlmostEqual(sh(components, -45, -90),  0)
		self.assertAlmostEqual(sh(components, -45, +90),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.86391580942710411)

	def test_sh_2_0(self) :
		components = self._fumaNormalizedComponents(2, 0)

		# axis
		self.assertAlmostEqual(sh(components,   0,   0), -.5)
		self.assertAlmostEqual(sh(components,   0, 180), -.5)
		self.assertAlmostEqual(sh(components,   0, -90), -.5)
		self.assertAlmostEqual(sh(components,   0, +90), -.5)
		self.assertAlmostEqual(sh(components, -90,   0),  +1)
		self.assertAlmostEqual(sh(components, +90,   0),  +1)
		# extremes
		self.assertAlmostEqual(sh(components,   0,  45), -.5)
		self.assertAlmostEqual(sh(components,   0, 135), -.5)
		self.assertAlmostEqual(sh(components,   0, 225), -.5)
		self.assertAlmostEqual(sh(components,   0, 315), -.5)

		self.assertAlmostEqual(sh(components,  45, -90), +.25)
		self.assertAlmostEqual(sh(components,  45, +90), +.25)
		self.assertAlmostEqual(sh(components, -45, -90), +.25)
		self.assertAlmostEqual(sh(components, -45, +90), +.25)

		self.assertAlmostEqual(sh(components,  45,   0), +.25)
		self.assertAlmostEqual(sh(components,  45, 180), +.25)
		self.assertAlmostEqual(sh(components, -45,   0), +.25)
		self.assertAlmostEqual(sh(components, -45, 180), +.25)

		self.assertAlmostEqual(sh(components,  30,   4),  -0.1250000000000011)

	def test_sh_2_m1(self) :
		components = self._fumaNormalizedComponents(2,-1)

		# axis
		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)
		self.assertAlmostEqual(sh(components,   0, -90),  0)
		self.assertAlmostEqual(sh(components,   0, +90),  0)
		# extremes
		self.assertAlmostEqual(sh(components,   0,  45),  0)
		self.assertAlmostEqual(sh(components,   0, 135),  0)
		self.assertAlmostEqual(sh(components,   0, 225),  0)
		self.assertAlmostEqual(sh(components,   0, 315),  0)

		self.assertAlmostEqual(sh(components,  45, -90), -1)
		self.assertAlmostEqual(sh(components,  45, +90), +1)
		self.assertAlmostEqual(sh(components, -45, -90), +1)
		self.assertAlmostEqual(sh(components, -45, +90), -1)

		self.assertAlmostEqual(sh(components,  45,   0),  0)
		self.assertAlmostEqual(sh(components,  45, 180),  0)
		self.assertAlmostEqual(sh(components, -45,   0),  0)
		self.assertAlmostEqual(sh(components, -45, 180),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.060410878340834695)


	def test_sh_2_m2(self) :
		components = self._fumaNormalizedComponents(2,-2)

		# axis
		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)
		self.assertAlmostEqual(sh(components,   0, -90),  0)
		self.assertAlmostEqual(sh(components,   0, +90),  0)
		# extremes
		self.assertAlmostEqual(sh(components,   0,  45), +1)
		self.assertAlmostEqual(sh(components,   0, 135), -1)
		self.assertAlmostEqual(sh(components,   0, 225), +1)
		self.assertAlmostEqual(sh(components,   0, 315), -1)

		self.assertAlmostEqual(sh(components,  45, -90),  0)
		self.assertAlmostEqual(sh(components,  45, +90),  0)
		self.assertAlmostEqual(sh(components, -45, -90),  0)
		self.assertAlmostEqual(sh(components, -45, +90),  0)

		self.assertAlmostEqual(sh(components,  45,   0),  0)
		self.assertAlmostEqual(sh(components,  45, 180),  0)
		self.assertAlmostEqual(sh(components, -45,   0),  0)
		self.assertAlmostEqual(sh(components, -45, 180),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.10437982572004907)

	def test_sh_3_3(self) :
		components = self._fumaNormalizedComponents(3,+3)

		# axis
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)
		self.assertAlmostEqual(sh(components,   0,   0), +1)
		self.assertAlmostEqual(sh(components,   0,  30),  0)
		self.assertAlmostEqual(sh(components,   0,  60), -1)
		self.assertAlmostEqual(sh(components,   0,  90),  0)
		self.assertAlmostEqual(sh(components,   0, 120), +1)
		self.assertAlmostEqual(sh(components,   0, 150),  0)
		self.assertAlmostEqual(sh(components,   0, 180), -1)
		self.assertAlmostEqual(sh(components,   0, 210),  0)
		self.assertAlmostEqual(sh(components,   0, 240), +1)
		self.assertAlmostEqual(sh(components,   0, 270),  0)
		self.assertAlmostEqual(sh(components,   0, 300), -1)
		self.assertAlmostEqual(sh(components,   0, 330),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.63532550316470549)

	def test_sh_3_2(self) :
		components = self._fumaNormalizedComponents(3,+2)

		# axis
		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)
		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0,  30),  0)
		self.assertAlmostEqual(sh(components,   0,  60),  0)
		self.assertAlmostEqual(sh(components,   0,  90),  0)
		self.assertAlmostEqual(sh(components,   0, 120),  0)
		self.assertAlmostEqual(sh(components,   0, 150),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components,   0, 210),  0)
		self.assertAlmostEqual(sh(components,   0, 240),  0)
		self.assertAlmostEqual(sh(components,   0, 270),  0)
		self.assertAlmostEqual(sh(components,   0, 300),  0)
		self.assertAlmostEqual(sh(components,   0, 330),  0)

		maxangle = math.degrees(math.asin(math.sqrt(1./3))) # 35.2644

		self.assertAlmostEqual(sh(components, +maxangle,   0), +1)
		self.assertAlmostEqual(sh(components, +maxangle,  45),  0)
		self.assertAlmostEqual(sh(components, +maxangle,  90), -1)
		self.assertAlmostEqual(sh(components, +maxangle, 135),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 180), +1)
		self.assertAlmostEqual(sh(components, +maxangle, 225),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 270), -1)
		self.assertAlmostEqual(sh(components, +maxangle, 315),  0)

		self.assertAlmostEqual(sh(components, -maxangle,   0), -1)
		self.assertAlmostEqual(sh(components, -maxangle,  45),  0)
		self.assertAlmostEqual(sh(components, -maxangle,  90), +1)
		self.assertAlmostEqual(sh(components, -maxangle, 135),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 180), -1)
		self.assertAlmostEqual(sh(components, -maxangle, 225),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 270), +1)
		self.assertAlmostEqual(sh(components, -maxangle, 315),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.96479696709759877)

	def test_sh_3_1(self) :
		components = self._fumaNormalizedComponents(3,+1)

		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)

		xymaxvalue = math.sqrt(3*5)*3/16 # -0.72618437741389052

		self.assertAlmostEqual(sh(components,   0,   0), -xymaxvalue)
		self.assertAlmostEqual(sh(components,   0,  90),  0)
		self.assertAlmostEqual(sh(components,   0, 180), +xymaxvalue)
		self.assertAlmostEqual(sh(components,   0, 270),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(1./5)))

		self.assertAlmostEqual(sh(components, +zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, +zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 270),  0)

		self.assertAlmostEqual(sh(components, -zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, -zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 270),  0)

		maxangle = math.degrees(math.acos(math.sqrt(4./15))) # 58.90907

		self.assertAlmostEqual(sh(components, +maxangle,   0), +1)
		self.assertAlmostEqual(sh(components, +maxangle,  90),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 180), -1)
		self.assertAlmostEqual(sh(components, +maxangle, 270),  0)

		self.assertAlmostEqual(sh(components, -maxangle,   0), +1)
		self.assertAlmostEqual(sh(components, -maxangle,  90),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 180), -1)
		self.assertAlmostEqual(sh(components, -maxangle, 270),  0)

		self.assertAlmostEqual(sh(components,  30,   4),  0.15684054105170947)

	def test_sh_3_0(self) :
		components = self._fumaNormalizedComponents(3, 0)

		# axis
		self.assertAlmostEqual(sh(components, -90,   0), -1)
		self.assertAlmostEqual(sh(components, +90,   0), +1)
		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0,  30),  0)
		self.assertAlmostEqual(sh(components,   0,  60),  0)
		self.assertAlmostEqual(sh(components,   0,  90),  0)
		self.assertAlmostEqual(sh(components,   0, 120),  0)
		self.assertAlmostEqual(sh(components,   0, 150),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components,   0, 210),  0)
		self.assertAlmostEqual(sh(components,   0, 240),  0)
		self.assertAlmostEqual(sh(components,   0, 270),  0)
		self.assertAlmostEqual(sh(components,   0, 300),  0)
		self.assertAlmostEqual(sh(components,   0, 330),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(3./5))) # 50.7684795164

		self.assertAlmostEqual(sh(components, +zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, +zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 270),  0)

		self.assertAlmostEqual(sh(components, -zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, -zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 270),  0)

		maxangle = math.degrees(math.asin(math.sqrt(3./15))) # 26.5650511771
		maxvalue = math.sqrt(3./15)

		self.assertAlmostEqual(sh(components, +maxangle,   0), -maxvalue)
		self.assertAlmostEqual(sh(components, +maxangle,  90), -maxvalue)
		self.assertAlmostEqual(sh(components, +maxangle, 180), -maxvalue)
		self.assertAlmostEqual(sh(components, +maxangle, 270), -maxvalue)

		self.assertAlmostEqual(sh(components, -maxangle,   0), +maxvalue)
		self.assertAlmostEqual(sh(components, -maxangle,  90), +maxvalue)
		self.assertAlmostEqual(sh(components, -maxangle, 180), +maxvalue)
		self.assertAlmostEqual(sh(components, -maxangle, 270), +maxvalue)

		self.assertAlmostEqual(sh(components,  30,   4),  -0.43750000000000006)

	def test_sh_3_m1(self) :
		components = self._fumaNormalizedComponents(3,-1)

		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)

		xymaxvalue = math.sqrt(3*5)*3/16 # -0.72618437741389052
		self.assertAlmostEqual(sh(components,   0,  90), -xymaxvalue)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components,   0, 270), +xymaxvalue)
		self.assertAlmostEqual(sh(components,   0,   0),  0)

		zeroangle = math.degrees(math.asin(math.sqrt(1./5)))

		self.assertAlmostEqual(sh(components, +zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, +zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, +zeroangle, 270),  0)

		self.assertAlmostEqual(sh(components, -zeroangle,   0),  0)
		self.assertAlmostEqual(sh(components, -zeroangle,  90),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 180),  0)
		self.assertAlmostEqual(sh(components, -zeroangle, 270),  0)

		maxangle = math.degrees(math.acos(math.sqrt(4./15))) # 58.90907
		self.assertAlmostEqual(sh(components, +maxangle,   0),  0)
		self.assertAlmostEqual(sh(components, +maxangle,  90), +1)
		self.assertAlmostEqual(sh(components, +maxangle, 180),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 270), -1)

		self.assertAlmostEqual(sh(components, -maxangle,   0),  0)
		self.assertAlmostEqual(sh(components, -maxangle,  90), +1)
		self.assertAlmostEqual(sh(components, -maxangle, 180),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 270), -1)

		self.assertAlmostEqual(sh(components,  30,   4),  0.010967359019241315)

	def test_sh_3_m2(self) :
		components = self._fumaNormalizedComponents(3,-2)

		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)

		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0,  30),  0)
		self.assertAlmostEqual(sh(components,   0,  60),  0)
		self.assertAlmostEqual(sh(components,   0,  90),  0)
		self.assertAlmostEqual(sh(components,   0, 120),  0)
		self.assertAlmostEqual(sh(components,   0, 150),  0)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components,   0, 210),  0)
		self.assertAlmostEqual(sh(components,   0, 240),  0)
		self.assertAlmostEqual(sh(components,   0, 270),  0)
		self.assertAlmostEqual(sh(components,   0, 300),  0)
		self.assertAlmostEqual(sh(components,   0, 330),  0)


		maxangle = math.degrees(math.asin(math.sqrt(1./3))) # 35.2644

		self.assertAlmostEqual(sh(components, +maxangle,   0),  0)
		self.assertAlmostEqual(sh(components, +maxangle,  45), +1)
		self.assertAlmostEqual(sh(components, +maxangle,  90),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 135), -1)
		self.assertAlmostEqual(sh(components, +maxangle, 180),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 225), +1)
		self.assertAlmostEqual(sh(components, +maxangle, 270),  0)
		self.assertAlmostEqual(sh(components, +maxangle, 315), -1)

		self.assertAlmostEqual(sh(components, -maxangle,   0),  0)
		self.assertAlmostEqual(sh(components, -maxangle,  45), -1)
		self.assertAlmostEqual(sh(components, -maxangle,  90),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 135), +1)
		self.assertAlmostEqual(sh(components, -maxangle, 180),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 225), -1)
		self.assertAlmostEqual(sh(components, -maxangle, 270),  0)
		self.assertAlmostEqual(sh(components, -maxangle, 315), +1)

		self.assertAlmostEqual(sh(components,  30,   4),  0.13559337107423225)

	def test_sh_3_m3(self) :
		components = self._fumaNormalizedComponents(3,-3)

		self.assertAlmostEqual(sh(components, -90,   0),  0)
		self.assertAlmostEqual(sh(components, +90,   0),  0)

		self.assertAlmostEqual(sh(components,   0,   0),  0)
		self.assertAlmostEqual(sh(components,   0,  30), +1)
		self.assertAlmostEqual(sh(components,   0,  60),  0)
		self.assertAlmostEqual(sh(components,   0,  90), -1)
		self.assertAlmostEqual(sh(components,   0, 120),  0)
		self.assertAlmostEqual(sh(components,   0, 150), +1)
		self.assertAlmostEqual(sh(components,   0, 180),  0)
		self.assertAlmostEqual(sh(components,   0, 210), -1)
		self.assertAlmostEqual(sh(components,   0, 240),  0)
		self.assertAlmostEqual(sh(components,   0, 270), +1)
		self.assertAlmostEqual(sh(components,   0, 300),  0)
		self.assertAlmostEqual(sh(components,   0, 330), -1)

		self.assertAlmostEqual(sh(components,  30,   4),  0.13504260449396654)

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

