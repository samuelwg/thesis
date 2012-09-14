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

"""
TODO:
- Peak/Zero testing of 4th order
+ Peak/Zero testing of 5th order
- Furse-Malham factors of 4th order
- Furse-Malham factors of 5th order
- SH Expansion of a function
	- 
- Resynthesis of a function
"""


from sphericalHarmonics import *
import unittest


class SHComponentIndexingAndEnumerationTests(unittest.TestCase) :

	def test_shIndexes_enumeratesIndexesUpToOrderN(self) :
		self.assertEqual(
			[
				(0,0),
				(1,-1), (1,0), (1,1),
				(2,-2), (2,-1), (2,0), (2,1), (2,2)
			],
			shIndexes(2),
			)

	def test_shi_mapsSphericalHarmonicsIndex2SquareMatrix(self) :

		self.assertEqual((0,0), shi(0, 0))

		self.assertEqual((1,0), shi(1,-1))
		self.assertEqual((1,1), shi(1, 0))
		self.assertEqual((0,1), shi(1,+1))

		self.assertEqual((2,0), shi(2,-2))
		self.assertEqual((2,1), shi(2,-1))
		self.assertEqual((2,2), shi(2, 0))
		self.assertEqual((1,2), shi(2,+1))
		self.assertEqual((0,2), shi(2,+2))


	def test_shi_reverse_mapsSquareMatrixIndex2SphericalHarmonics(self) :

		self.assertEqual(shi_reverse(0,0), (0, 0))

		self.assertEqual(shi_reverse(1,0), (1,-1))
		self.assertEqual(shi_reverse(1,1), (1, 0))
		self.assertEqual(shi_reverse(0,1), (1,+1))

		self.assertEqual(shi_reverse(2,0), (2,-2))
		self.assertEqual(shi_reverse(2,1), (2,-1))
		self.assertEqual(shi_reverse(2,2), (2, 0))
		self.assertEqual(shi_reverse(1,2), (2,+1))
		self.assertEqual(shi_reverse(0,2), (2,+2))


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


"""
Peaks obtained from the following code:

	import sympy as sp
	from sympy.abc import z
	for l in xrange(5+1) :
		for m in xrange(l+1) :
			try:
				print "shExtremes(%i,%i) ="%(l,m), solve(
					simplify(diff(assoc_legendre(l,m,z),z)), z)
			except NotImplementedError :
				print "Unsolvable", sp.diff(sp.assoc_legendre(l,m,z),z)

For some reason 4,1 is not solveable so, it has been obtained by hand.
Also the former code does not include z=+1,-1 for most l,0 components
because these extremes have no zero derivative, still growing but the
slope is not zero yet.
"""
from sympy import sqrt, S, simplify
shExtremes = [
((0,0) , []),
((1,0) , []), # 1,-1
((1,1) , [0]),
((2,0) , [0]), # 1,-1
((2,1) , [1/sqrt(2)]),
((2,2) , [0]),
((3,0) , [ # 1,-1
			+1/sqrt(5),
			-1/sqrt(5),
		]),
((3,1) , [
			0,
			-sqrt(11)/sqrt(15),
			+sqrt(11)/sqrt(15),
		]),
((3,2) , [
			+1/sqrt(3),
			-1/sqrt(3),
		]),
((3,3) , [ # +1 -1 are false extrema
			0,
			+1,
			-1,
		]),
((4,0) , [ # 1,-1
			-sqrt(3)/sqrt(7),
			+sqrt(3)/sqrt(7),
			0,
		]),
((4,1) , []), # Unsolvable by numpy
((4,2) , [
			+2/sqrt(7),
			0,
			-2/sqrt(7),
		]),
((4,3) , [
			+1,
			-1,
			-S(1)/2,
			+S(1)/2,
			]),
((4,4) , [ # +1 -1 are false extrema
			+1,
			-1,
			0,
		]),
((5,0) , [ # 1,-1
			+sqrt(S(1)/3 +2*sqrt(7)/21 ),
			-sqrt(S(1)/3 +2*sqrt(7)/21 ),
			+sqrt(S(1)/3 -2*sqrt(7)/21 ),
			-sqrt(S(1)/3 -2*sqrt(7)/21 ),
		]),
((5,1) , [
			0,
			+sqrt(+63 +2*sqrt(231))/sqrt(105),
			-sqrt(+63 +2*sqrt(231))/sqrt(105),
			+sqrt(+63 -2*sqrt(231))/sqrt(105),
			-sqrt(+63 -2*sqrt(231))/sqrt(105),
		]),
((5,2) , [
			+sqrt(S(2)/5 +sqrt(S(7)/3)/5 ),
			-sqrt(S(2)/5 -sqrt(S(7)/3)/5 ),
			-sqrt(S(2)/5 +sqrt(S(7)/3)/5 ),
			+sqrt(S(2)/5 -sqrt(S(7)/3)/5 ),
		]),
((5,3) , [
			0,
			+1,
			-1,
			-sqrt(7)/sqrt(15),
			+sqrt(7)/sqrt(15),
		]),
((5,4) , [
			+1,
			+1/sqrt(5),
			-1,
			-1/sqrt(5),
		]),
((5,5) , [ # +1 -1 are false extrema
			0,
			+1,
			-1,
		]),
]

shExtremes2 = dict([
((0,0) , []),
((1,0) , []), # 1,-1
((1,1) , [0]),
((2,0) , [0]), # 1,-1
((2,1) , [1/sqrt(2)]),
((2,2) , [0]),
((3,0) , [ # 1,-1
			+1/sqrt(5),
			-1/sqrt(5),
		]),
((3,1) , [
			0,
			-sqrt(11)/sqrt(15),
			+sqrt(11)/sqrt(15),
		]),
((3,2) , [
			1/sqrt(3),
			-1/sqrt(3),
		]),
((3,3) , [ # +1 -1 are false extrema
			0,
			+1,
			-1,
		]),
((4,0) , [ # 1,-1
			-sqrt(3)/sqrt(7),
			+sqrt(3)/sqrt(7),
			0,
		]),
((4,1) , []), # Unsolvable by numpy
((4,2) , [
			+2/sqrt(7),
			0,
			-2/sqrt(7),
		]),
((4,3) , [
			+1,
			-1,
			-S(1)/2,
			+S(1)/2,
			]),
((4,4) , [ # +1 -1 are false extrema
			+1,
			-1,
			0,
		]),
((5,0) , [ # 1,-1
			+sqrt((1 +2/sqrt(7))/3),
			-sqrt((1 +2/sqrt(7))/3),
			+sqrt((1 -2/sqrt(7))/3),
			-sqrt((1 -2/sqrt(7))/3),
		]),
((5,1) , [
			0,
			+sqrt(+3 +2*sqrt(11)/sqrt(3)/sqrt(7))/sqrt(5),
			-sqrt(+3 +2*sqrt(11)/sqrt(3)/sqrt(7))/sqrt(5),
			+sqrt(+3 -2*sqrt(11)/sqrt(3)/sqrt(7))/sqrt(5),
			-sqrt(+3 -2*sqrt(11)/sqrt(3)/sqrt(7))/sqrt(5),
		]),
((5,2) , [
			+sqrt((2 +sqrt(S(7)/3))/5),
			-sqrt((2 -sqrt(S(7)/3))/5),
			-sqrt((2 +sqrt(S(7)/3))/5),
			+sqrt((2 -sqrt(S(7)/3))/5),
		]),
((5,3) , [
			0,
			+1,
			-1,
			-sqrt(7)/sqrt(15),
			+sqrt(7)/sqrt(15),
		]),
((5,4) , [
			+1,
			+1/sqrt(5),
			-1,
			-1/sqrt(5),
		]),
((5,5) , [ # +1 -1 are false extrema
			0,
			+1,
			-1,
		]),
])

shZerosRef = [
((0,0) , []),
((1,0) , [0]),
((1,1) , [1, -1]),
((2,0) , [3**(S(1)/2)/3, -3**(S(1)/2)/3]),
((2,1) , [0, 1, -1]),
((2,2) , [1, -1]),
((3,0) , [-15**(S(1)/2)/5, 15**(S(1)/2)/5, 0]),
((3,1) , [1, -1, 5**(S(1)/2)/5, -5**(S(1)/2)/5]),
((3,2) , [1, -1, 0]),
((3,3) , [1, -1, 1, -1]),
((4,0) , [(2*30**(S(1)/2)/35 + S(3)/7)**(S(1)/2), -(-2*30**(S(1)/2)/35 + S(3)/7)**(S(1)/2), -(2*30**(S(1)/2)/35 + S(3)/7)**(S(1)/2), (-2*30**(S(1)/2)/35 + S(3)/7)**(S(1)/2)]),
((4,1) , [0, 1, -1, -21**(S(1)/2)/7, 21**(S(1)/2)/7]),
((4,2) , [1, 7**(S(1)/2)/7, -1, -7**(S(1)/2)/7]),
((4,3) , [0, 1, -1, 1, -1]),
((4,4) , [1, -1]),
((5,0) , [-(2*70**(S(1)/2)/63 + S(5)/9)**(S(1)/2), (2*70**(S(1)/2)/63 + S(5)/9)**(S(1)/2), -(-2*70**(S(1)/2)/63 + S(5)/9)**(S(1)/2), 0, (-2*70**(S(1)/2)/63 + S(5)/9)**(S(1)/2)]),
((5,1) , [1, -1, 21**(S(1)/2)*(2*7**(S(1)/2) + 7)**(S(1)/2)/21, -21**(S(1)/2)*(2*7**(S(1)/2) + 7)**(S(1)/2)/21, 21**(S(1)/2)*(-2*7**(S(1)/2) + 7)**(S(1)/2)/21, -21**(S(1)/2)*(-2*7**(S(1)/2) + 7)**(S(1)/2)/21]),
((5,2) , [1, 3**(S(1)/2)/3, -1, 0, -3**(S(1)/2)/3]),
((5,3) , [1, -1, 1, S(1)/3, -1, -S(1)/3]),
((5,4) , [1, -1, 0]),
((5,5) , [1, -1]),
]


shZeros = dict([
((0,0) , []),
((1,0) , [0]),
((1,1) , [1, -1]),
((2,0) , [
			+1/sqrt(3),
			-1/sqrt(3),
		]),
((2,1) , [0, 1, -1]),
((2,2) , [1, -1]),
((3,0) , [
			-sqrt(S(3)/5),
			+sqrt(S(3)/5),
			0,
		]),
((3,1) , [
			1,
			-1,
			5**(S(1)/2)/5,
			-5**(S(1)/2)/5],
			),
((3,2) , [1, -1, 0]),
((3,3) , [1, -1, 1, -1]),
((4,0) , [
			+sqrt((+2*sqrt(S(6)/5) + S(3))/7),
			-sqrt((-2*sqrt(S(6)/5) + S(3))/7),
			-sqrt((+2*sqrt(S(6)/5) + S(3))/7),
			+sqrt((-2*sqrt(S(6)/5) + S(3))/7),
		]),
((4,1) , [
			0,
			1,
			-1,
			-sqrt(S(3)/7),
			+sqrt(S(3)/7),
		]),
((4,2) , [
			+1,
			+1/sqrt(7),
			-1,
			-1/sqrt(7),
		]),
((4,3) , [0, 1, -1, 1, -1]),
((4,4) , [1, -1]),
((5,0) , [
			-sqrt(S(5) +2*sqrt(S(10)/7))/3,
			+sqrt(S(5) +2*sqrt(S(10)/7))/3,
			-sqrt(S(5) -2*sqrt(S(10)/7))/3,
			0,
			+sqrt(S(5) -2*sqrt(S(10)/7))/3,
		]),
((5,1) , [
			+1,
			-1,
			+sqrt(S(1)/3 +S(2)/3/sqrt(7)),
			-sqrt(S(1)/3 +S(2)/3/sqrt(7)),
			+sqrt(S(1)/3 -S(2)/3/sqrt(7)),
			-sqrt(S(1)/3 -S(2)/3/sqrt(7)),
		]),
((5,2) , [
			1,
			+1/sqrt(3),
			-1,
			0,
			-1/sqrt(3),
		]),
((5,3) , [1, -1, 1, S(1)/3, -1, -S(1)/3]),
((5,4) , [1, -1, 0]),
((5,5) , [1, -1]),
])


class NodalSHStructureTests(unittest.TestCase) :
	def testMaxMins(self) :
		for ((l,m), zs) in shExtremes :
			otherZs = shExtremes2[l,m]
			self.assertEqual(len(zs), len(otherZs),
				"l=%i,m=%i differ in length, expected %i, result %i"%(
					l,m,len(zs), len(otherZs)))
			for z1,z2 in zip(zs, otherZs) :
				self.assertEqual(simplify(z1 - z2), 0)

	def testZeros(self) :
		for ((l,m), zs) in shZerosRef :
			otherZs = shZeros[l,m]
			self.assertEqual(len(zs), len(otherZs),
				"l=%i,m=%i differ in length, expected %i, result %i"%(
					l,m,len(zs), len(otherZs)))
			for z1,z2 in zip(zs, otherZs) :
				self.assertEqual(simplify(z1 - z2), 0)





class SphericalHarmonicsTests(unittest.TestCase) :

	def sh(self, e, a) :
		"""
		Decodes the sh components for (e,a) position.
		"""
		return (self._components*semiNormalizedSH(e,a)).sum()

	def assertEqualAt(self, ea, value) :
		e,a=ea
		self.assertAlmostEqual(self.sh(e, a), value)

	def prepareOrder(self, l,m) :
		# Create projection of a single SH component,
		# normalized alla Furse-Malham so that the result
		# max values of 1. in max but for 0,0
		self._components = np.zeros(shShape)
		self._components[shi(l,m)] = fuma[shi(l,m)]



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

		from math import asin, sqrt, degrees

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

		maxangle = degrees(asin(sqrt(1./3))) # 35.2644

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

		from math import asin, sqrt, degrees

		self.assertEqualAt(( -90,   0),  0)
		self.assertEqualAt(( +90,   0),  0)

		xymaxvalue = sqrt(3*5)*3/16 # -0.72618437741389052

		self.assertEqualAt((   0,   0), -xymaxvalue)
		self.assertEqualAt((   0,  90),  0)
		self.assertEqualAt((   0, 180), +xymaxvalue)
		self.assertEqualAt((   0, 270),  0)

		zeroangle = degrees(asin(sqrt(1./5)))

		self.assertEqualAt(( +zeroangle,   0),  0)
		self.assertEqualAt(( +zeroangle,  90),  0)
		self.assertEqualAt(( +zeroangle, 180),  0)
		self.assertEqualAt(( +zeroangle, 270),  0)

		self.assertEqualAt(( -zeroangle,   0),  0)
		self.assertEqualAt(( -zeroangle,  90),  0)
		self.assertEqualAt(( -zeroangle, 180),  0)
		self.assertEqualAt(( -zeroangle, 270),  0)

		maxangle = degrees(asin(sqrt(11./15))) # 58.90907

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

		from math import asin, sqrt, degrees

		zeroangle = degrees(asin(sqrt(3./5))) # 50.7684795164
		maxangle = degrees(asin(sqrt(1./5))) # 26.5650511771
		maxvalue = sqrt(3./15)

		self.assertPeakAtElevation( 0, +90, +1)
		self.assertPeakAtElevation( 0, +zeroangle, 0)
		self.assertPeakAtElevation( 0, +maxangle, -maxvalue)
		self.assertPeakAtElevation( 0, 0, 0)
		self.assertPeakAtElevation( 0, -maxangle, +maxvalue)
		self.assertPeakAtElevation( 0, -zeroangle, 0)
		self.assertPeakAtElevation( 0, -90, -1)

		self.assertEqualAt((  30,   4),  -0.43750000000000006)

	def test_sh_3_m1(self) :
		self.prepareOrder(3,-1)

		from math import asin, sqrt, degrees

		xymaxvalue = sqrt(3*5)*3/16 # -0.72618437741389052
		zeroangle = degrees(asin(sqrt(1./5)))
		maxangle = degrees(asin(sqrt(11./15))) # 58.90907

		self.assertPeakAtElevation(-1, +90, 0)
		self.assertPeakAtElevation(-1, +maxangle, +1)
		self.assertPeakAtElevation(-1, +zeroangle, 0)
		self.assertPeakAtElevation(-1, 0, -xymaxvalue)
		self.assertPeakAtElevation(-1, -zeroangle, 0)
		self.assertPeakAtElevation(-1, -maxangle, +1)
		self.assertPeakAtElevation(-1, -90, 0)

		self.assertEqualAt((  30,   4),  0.010967359019241315)

	def test_sh_3_m2(self) :
		self.prepareOrder(3,-2)

		from math import asin, sqrt, degrees

		maxangle = degrees(asin(sqrt(1./3))) # 35.2644

		self.assertPeakAtElevation(-2, +90, 0)
		self.assertPeakAtElevation(-2, +maxangle, +1)
		self.assertPeakAtElevation(-2,   0, 0)
		self.assertPeakAtElevation(-2, -maxangle, -1)
		self.assertPeakAtElevation(-2, -90, 0)

		self.assertEqualAt((  30,   4),  0.13559337107423225)

	def test_sh_3_m3(self) :
		self.prepareOrder(3,-3)

		self.assertPeakAtElevation(-3, +90, 0)
		self.assertPeakAtElevation(-3,   0, 1)
		self.assertPeakAtElevation(-3, -90, 0)

		self.assertEqualAt((  30,   4),  0.13504260449396654)

	# TODO: Warning 4th and above are not FuMa normalized yet!!!

	def test_sh_4_0(self) :
		self.prepareOrder(4,-0)

		from math import asin, sqrt, degrees

		zero2     = degrees(asin( sqrt(3./7 + 2*sqrt(6./5)/7) ))
		maxangle1 = degrees(asin( sqrt(3./7) )) # 40.8933946491309
		zero1     = degrees(asin( sqrt(3./7 - 2*sqrt(6./5)/7) ))
		maxvalue1 = 3./7
		maxvalue0 = 3./8

		self.assertPeakAtElevation( 0, +90, +1)
		self.assertPeakAtElevation( 0, +zero2, 0)
		self.assertPeakAtElevation( 0, +maxangle1, -maxvalue1)
		self.assertPeakAtElevation( 0, +zero1, 0)
		self.assertPeakAtElevation( 0,  0, maxvalue0)
		self.assertPeakAtElevation( 0, -zero1, 0)
		self.assertPeakAtElevation( 0, -maxangle1, -maxvalue1)
		self.assertPeakAtElevation( 0, -zero2, 0)
		self.assertPeakAtElevation( 0, -90, +1)

		self.assertEqualAt((  30,   4), -0.28906249999999989)

	def test_sh_4_1(self) :
		self.prepareOrder(4,+1)

		from math import asin, sqrt, degrees

		zeroangle = degrees(asin(sqrt(3./7)))
		maxangle1 = degrees(asin(sqrt( (27-sqrt(393))/7/8) ))
		maxvalue1 = 0.55571119769376354 # b2b, not analytically checked
		maxangle2 = degrees(asin(sqrt( (27+sqrt(393))/7/8) )) # 66.1221762567095
		maxvalue2 = 0.83486203359622035 # b2b, not analytically checked

		self.assertPeakAtElevation(+1, +90, 0)
		self.assertPeakAtElevation(+1, +maxangle2, +maxvalue2)
		self.assertPeakAtElevation(+1, +zeroangle, 0)
		self.assertPeakAtElevation(+1, +maxangle1, -maxvalue1)
		self.assertPeakAtElevation(+1,  0, 0)
		self.assertPeakAtElevation(+1, -maxangle1, +maxvalue1)
		self.assertPeakAtElevation(+1, -zeroangle, 0)
		self.assertPeakAtElevation(+1, -maxangle2, -maxvalue2)
		self.assertPeakAtElevation(+1, -90, 0)

		self.assertEqualAt((  30,   4), -0.42686588506525242)

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

	def assertPeakAtElevation(self, degree, elevation, maxvalue) :
		if degree > 0 :
			for i in xrange(degree) :
				self.assertEqualAt((elevation, 90.*(i*4+0)/degree), +maxvalue)
				self.assertEqualAt((elevation, 90.*(i*4+1)/degree),  0)
				self.assertEqualAt((elevation, 90.*(i*4+2)/degree), -maxvalue)
				self.assertEqualAt((elevation, 90.*(i*4+3)/degree),  0)
		elif degree < 0 :
			degree = abs(degree)
			for i in xrange(degree) :
				self.assertEqualAt((elevation, 90.*(i*4+0)/degree),  0)
				self.assertEqualAt((elevation, 90.*(i*4+1)/degree), +maxvalue)
				self.assertEqualAt((elevation, 90.*(i*4+2)/degree),  0)
				self.assertEqualAt((elevation, 90.*(i*4+3)/degree), -maxvalue)
		else :
			subdivision = 19
			for i in xrange(subdivision) :
				self.assertEqualAt((elevation, 360.*(i)/subdivision),  maxvalue)

	def test_sh_5_m5(self) :
		self.prepareOrder(5,-5)

		from math import sqrt, asin, degrees

		maxvalue = 3*sqrt(7./2)/8 # 0.70156076002011403

		self.assertPeakAtElevation(-5, +90, 0)
		self.assertPeakAtElevation(-5, 0, maxvalue)
		self.assertPeakAtElevation(-5, -90, 0)

		self.assertEqualAt((  30,   4), +0.116888055250392)


	def test_sh_5_m4(self) :
		self.prepareOrder(5,-4)

		from math import sqrt, asin, degrees

		peakelevation = degrees(asin( 1/sqrt(5) ))
		maxvalue = 6./25*sqrt(7)

		self.assertPeakAtElevation(-4, +90, 0)
		self.assertPeakAtElevation(-4, +peakelevation, +maxvalue)
		self.assertPeakAtElevation(-4,   0, 0)
		self.assertPeakAtElevation(-4, -peakelevation, -maxvalue)
		self.assertPeakAtElevation(-4, -90, 0)

		self.assertEqualAt((  30,   4), +0.171987108913405)

	def test_sh_5_m3(self) :
		self.prepareOrder(5,-3)

		from math import sqrt, asin, degrees

		peakelevation  = degrees(asin( sqrt(7./15) ))
		zeroelevation = degrees(asin( 1./3 ))
		maxvalue = 0.65174409883816409
		maxvalueAtZero = 0.52291251658379723

		self.assertPeakAtElevation(-3, +90, 0)
		self.assertPeakAtElevation(-3, +peakelevation, +maxvalue)
		self.assertPeakAtElevation(-3, +zeroelevation, 0)
		self.assertPeakAtElevation(-3, 0, -maxvalueAtZero)
		self.assertPeakAtElevation(-3, -zeroelevation, 0)
		self.assertPeakAtElevation(-3, -peakelevation, +maxvalue)
		self.assertPeakAtElevation(-3, -90, 0)

		self.assertEqualAt((  30,   4), +0.0882693352024631)

	def test_sh_5_m2(self) :
		self.prepareOrder(5,-2)

		from math import sqrt, asin, degrees

		peakelevation2 = degrees(asin( sqrt((2 +sqrt(7./3))/5) ))
		zeroelevation = degrees(asin( sqrt(1./3) ))
		peakelevation1 = degrees(asin( sqrt((2 -sqrt(7./3))/5) ))
		maxvalue1 = 0.51092271611124063
		maxvalue2 = 0.70750122131450299

		self.assertPeakAtElevation(-2, -90, 0)
		self.assertPeakAtElevation(-2, +peakelevation2, +maxvalue2)
		self.assertPeakAtElevation(-2, +zeroelevation, 0)
		self.assertPeakAtElevation(-2, +peakelevation1, -maxvalue1)
		self.assertPeakAtElevation(-2, 0, 0)
		self.assertPeakAtElevation(-2, -peakelevation1, +maxvalue1)
		self.assertPeakAtElevation(-2, -zeroelevation, 0)
		self.assertPeakAtElevation(-2, -peakelevation2, -maxvalue2)
		self.assertPeakAtElevation(-2, +90, 0)

		self.assertEqualAt((  30,   4), -0.0334242167222746)

	def test_sh_5_m1(self) :
		self.prepareOrder(5,-1)

		from math import sqrt, asin, degrees

		peakelevation2  = degrees(asin( sqrt( +3 +2*sqrt(11)/sqrt(3*7))/sqrt(5) ))
		peakelevation1  = degrees(asin( sqrt( +3 -2*sqrt(11)/sqrt(3*7))/sqrt(5) ))
		zeroelevation2 = degrees(asin( sqrt( (1+2/sqrt(7))/3) ))
		zeroelevation1 = degrees(asin( sqrt( (1-2/sqrt(7))/3) ))
		maxvalue2 = 0.83078704673985981
		maxvalue1 = 0.53159466040326442
		maxvalue0 = 0.48412291827592713

		self.assertPeakAtElevation(-1, -90, 0)
		self.assertPeakAtElevation(-1, +peakelevation2, +maxvalue2)
		self.assertPeakAtElevation(-1, +zeroelevation2, 0)
		self.assertPeakAtElevation(-1, +peakelevation1, -maxvalue1)
		self.assertPeakAtElevation(-1, +zeroelevation1, 0)
		self.assertPeakAtElevation(-1,               0, +maxvalue0)
		self.assertPeakAtElevation(-1, -zeroelevation1, 0)
		self.assertPeakAtElevation(-1, -peakelevation1, -maxvalue1)
		self.assertPeakAtElevation(-1, -zeroelevation2, 0)
		self.assertPeakAtElevation(-1, -peakelevation2, +maxvalue2)
		self.assertPeakAtElevation(-1, +90, 0)

		self.assertEqualAt((  30,   4), -0.0347299702275976)


	def test_sh_5_0(self) :
		self.prepareOrder(5, 0)

		from math import sqrt, degrees, asin
		zeroelevation2 = degrees(asin( sqrt(5 +2*sqrt(10./7))/3 ))
		zeroelevation1 = degrees(asin( sqrt(5 -2*sqrt(10./7))/3 ))
		peakelevation2 = degrees(asin( sqrt((1 +2/sqrt(7))/3) ))
		peakelevation1 = degrees(asin( sqrt((1 -2/sqrt(7))/3) ))
		peakelevation0  = 0
		maxvalue1 = 0.34662772505541772
		maxvalue2 = 0.41969693413128728

		self.assertPeakAtElevation(0, +90,+1)
		self.assertPeakAtElevation(0, +zeroelevation2, 0)
		self.assertPeakAtElevation(0, +peakelevation2, -maxvalue2)
		self.assertPeakAtElevation(0, +zeroelevation1, 0)
		self.assertPeakAtElevation(0, +peakelevation1, +maxvalue1)
		self.assertPeakAtElevation(0,               0, +0)
		self.assertPeakAtElevation(0, -peakelevation1, -maxvalue1)
		self.assertPeakAtElevation(0, -zeroelevation1, 0)
		self.assertPeakAtElevation(0, -peakelevation2 ,+maxvalue2)
		self.assertPeakAtElevation(0, -zeroelevation2, 0)
		self.assertPeakAtElevation(0, -90,-1)

		self.assertEqualAt((  30,   4), +0.0898437500000000)

	def test_sh_5_1(self) :
		self.prepareOrder(5,+1)

		from math import sqrt, asin, degrees

		peakelevation2  = degrees(asin( sqrt( +3 +2*sqrt(11)/sqrt(3*7))/sqrt(5) ))
		peakelevation1  = degrees(asin( sqrt( +3 -2*sqrt(11)/sqrt(3*7))/sqrt(5) ))
		zeroelevation2 = degrees(asin( sqrt( (1+2/sqrt(7))/3) ))
		zeroelevation1 = degrees(asin( sqrt( (1-2/sqrt(7))/3) ))
		maxvalue2 = 0.83078704673985981
		maxvalue1 = 0.53159466040326442
		maxvalue0 = 0.48412291827592713

		self.assertPeakAtElevation(+1, -90, 0)
		self.assertPeakAtElevation(+1, +peakelevation2 , +maxvalue2)
		self.assertPeakAtElevation(+1, +zeroelevation2, 0)
		self.assertPeakAtElevation(+1, +peakelevation1 , -maxvalue1)
		self.assertPeakAtElevation(+1, +zeroelevation1, 0)
		self.assertPeakAtElevation(+1,               0, +maxvalue0)
		self.assertPeakAtElevation(+1, -zeroelevation1, 0)
		self.assertPeakAtElevation(+1, -peakelevation1 , -maxvalue1)
		self.assertPeakAtElevation(+1, -zeroelevation2, 0)
		self.assertPeakAtElevation(+1, -peakelevation2 , +maxvalue2)
		self.assertPeakAtElevation(+1, +90, 0)

		self.assertEqualAt((  30,   4), -0.496661713330414)

	def test_sh_5_2(self) :
		self.prepareOrder(5,+2)

		from math import sqrt, asin, degrees

		peakelevation2 = degrees(asin( sqrt((2 +sqrt(7./3))/5) ))
		zeroelevation = degrees(asin( sqrt(1./3) ))
		peakelevation1 = degrees(asin( sqrt((2 -sqrt(7./3))/5) ))
		maxvalue1 = 0.51092271611124063
		maxvalue2 = 0.70750122131450299

		self.assertPeakAtElevation(+2, -90, 0)
		self.assertPeakAtElevation(+2, +peakelevation2, +maxvalue2)
		self.assertPeakAtElevation(+2, +zeroelevation, 0)
		self.assertPeakAtElevation(+2, +peakelevation1, -maxvalue1)
		self.assertPeakAtElevation(+2, 0, 0)
		self.assertPeakAtElevation(+2, -peakelevation1, +maxvalue1)
		self.assertPeakAtElevation(+2, -zeroelevation, 0)
		self.assertPeakAtElevation(+2, -peakelevation2, -maxvalue2)
		self.assertPeakAtElevation(+2, +90, 0)

		self.assertEqualAt((  30,   4), -0.237825659660080)

	def test_sh_5_3(self) :
		self.prepareOrder(5,+3)

		from math import sqrt, asin, degrees

		peakelevation  = degrees(asin( sqrt(7./15) ))
		zeroelevation = degrees(asin( 1./3 ))
		maxvalue = 0.65174409883816409
		maxvalueAtZero = 0.52291251658379723

		self.assertPeakAtElevation(+3, +90, 0)
		self.assertPeakAtElevation(+3, +peakelevation, +maxvalue)
		self.assertPeakAtElevation(+3, +zeroelevation, 0)
		self.assertPeakAtElevation(+3, 0, -maxvalueAtZero)
		self.assertPeakAtElevation(+3, -zeroelevation, 0)
		self.assertPeakAtElevation(+3, -peakelevation, +maxvalue)
		self.assertPeakAtElevation(+3, -90, 0)

		self.assertEqualAt((  30,   4), +0.415274572137154)

	def test_sh_5_4(self) :
		self.prepareOrder(5,+4)

		from math import sqrt, asin, degrees

		peakelevation = degrees(asin( 1/sqrt(5) ))
		maxvalue = 6./25*sqrt(7)

		self.assertPeakAtElevation(+4, +90, 0)
		self.assertPeakAtElevation(+4, +peakelevation, +maxvalue)
		self.assertPeakAtElevation(+4,   0, 0)
		self.assertPeakAtElevation(+4, -peakelevation, -maxvalue)
		self.assertPeakAtElevation(+4, -90, 0)

		self.assertEqualAt((  30,   4), +0.599790327779047)

	def test_sh_5_5(self) :
		self.prepareOrder(5,+5)

		from math import sqrt, asin, degrees

		maxvalue = 3*sqrt(7./2)/8 # 0.70156076002011403

		self.assertPeakAtElevation(+5, +90, 0)
		self.assertPeakAtElevation(+5, 0, maxvalue)
		self.assertPeakAtElevation(+5, -90, 0)

		self.assertEqualAt((  30,   4), +0.321147292404416)

@unittest.skip("Slow")
class SphericalHarmonicsTests_cartesianRecursive(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_cartesianRecursive(e,a)).sum()

@unittest.skip("Slow")
class SphericalHarmonicsTests_cartesian(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_cartesian(e,a)).sum()

@unittest.skip("Slow")
class SphericalHarmonicsTests_lambdified(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_lambdified(e,a)).sum()

@unittest.skip("Slow")
class SphericalHarmonicsTests_genericNumpy(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_genericNumpy(e,a)).sum()

@unittest.skip("Slow")
class SphericalHarmonicsTests_sympyGeneratedExpressions(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_sympyGeneratedExpressions(e,a)).sum()




if __name__ == "__main__" :
	sys.exit(unittest.main())

