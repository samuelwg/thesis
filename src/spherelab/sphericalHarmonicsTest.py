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
from sympy import sqrt, S
shExtremes = [
((0,0) , []),
((1,0) , []),
((1,1) , [0]),
((2,0) , [0]),
((2,1) , [sqrt(2)/2]),
((2,2) , [0]),
((3,0) , [sqrt(5)/5, -sqrt(5)/5]),
((3,1) , [0, -sqrt(165)/15, sqrt(165)/15]),
((3,2) , [sqrt(3)/3, -sqrt(3)/3]),
((3,3) , [0, 1, -1]),
((4,0) , [-sqrt(21)/7, sqrt(21)/7, 0]),
((4,1) , []), # Unsolvable
((4,2) , [2*sqrt(7)/7, 0, -2*sqrt(7)/7]),
((4,3) , [1, -1, -S(1)/2, S(1)/2]),
((4,4) , [1, -1, 0]),
((5,0) , [(2*sqrt(7)/21 + S(1)/3)**S(1./2), -(2*sqrt(7)/21 + S(1)/3)**S(1./2), (-2*sqrt(7)/21 + S(1)/3)**S(1./2), -(-2*sqrt(7)/21 + S(1)/3)**S(1./2)]),
((5,1) , [0, sqrt(105)*(2*sqrt(231) + 63)**S(1./2)/105, -sqrt(105)*(2*sqrt(231) + 63)**S(1./2)/105, sqrt(105)*(-2*sqrt(231) + 63)**S(1./2)/105, -sqrt(105)*(-2*sqrt(231) + 63)**S(1./2)/105]),
((5,2) , [(sqrt(21)/15 + S(2)/5)**S(1./2), -(-sqrt(21)/15 + S(2)/5)**S(1./2), -(sqrt(21)/15 + S(2)/5)**S(1./2), (-sqrt(21)/15 + S(2)/5)**S(1./2)]),
((5,3) , [0, 1, -1, -sqrt(105)/15, sqrt(105)/15]),
((5,4) , [1, sqrt(5)/5, -1, -sqrt(5)/5]),
((5,5) , [0, 1, -1]),
]

shExtremes2 = dict([
((0,0) , []),
((1,0) , []),
((1,1) , [0]),
((2,0) , [0]),
((2,1) , [1/sqrt(2)]),
((2,2) , [0]),
((3,0) , [1/sqrt(5), -1/sqrt(5)]),
((3,1) , [0, -sqrt(11)/sqrt(15), sqrt(S(11)/15)]),
((3,2) , [1/sqrt(3), -1/sqrt(3)]),
((3,3) , [0, 1, -1]),
((4,0) , [
			-sqrt(3)/sqrt(7),
			+sqrt(3)/sqrt(7),
			0,
		]),
((4,1) , []), # Unsolvable
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
((4,4) , [
			+1,
			-1,
			0,
		]),
((5,0) , [
			+(2*sqrt(7)/21 + S(1)/3)**S(1./2),
			-(2*sqrt(7)/21 + S(1)/3)**S(1./2),
			+(-2*sqrt(7)/21 + S(1)/3)**S(1./2),
			-(-2*sqrt(7)/21 + S(1)/3)**S(1./2)
		]),
((5,1) , [
			0,
			+sqrt(105)*(2*sqrt(231) + 63)**S(1./2)/105,
			-sqrt(105)*(2*sqrt(231) + 63)**S(1./2)/105,
			+sqrt(105)*(-2*sqrt(231) + 63)**S(1./2)/105,
			-sqrt(105)*(-2*sqrt(231) + 63)**S(1./2)/105,
		]),
((5,2) , [
			+sqrt(+sqrt(21)/15 + S(2)/5),
			-sqrt(-sqrt(21)/15 + S(2)/5),
			-sqrt(+sqrt(21)/15 + S(2)/5),
			+sqrt(-sqrt(21)/15 + S(2)/5),
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
((5,5) , [
			0,
			+1,
			-1,
		]),
])


class NodalSHStructureTests(unittest.TestCase) :
	def test(self) :
		""""""
		for ((l,m), zs) in shExtremes :
			self.assertEqual(zs, shExtremes2[l,m])





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

#@unittest.skip("Slow")
class SphericalHarmonicsTests_cartesianRecursive(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_cartesianRecursive(e,a)).sum()

#@unittest.skip("Slow")
class SphericalHarmonicsTests_cartesian(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_cartesian(e,a)).sum()

#@unittest.skip("Slow")
class SphericalHarmonicsTests_lambdified(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_lambdified(e,a)).sum()

#@unittest.skip("Slow")
class SphericalHarmonicsTests_genericNumpy(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_genericNumpy(e,a)).sum()

#@unittest.skip("Slow")
class SphericalHarmonicsTests_sympyGeneratedExpressions(SphericalHarmonicsTests) :
	def sh(self, e, a) :
		"""Decodes the sh components for (e,a) position."""
		return (self._components*semiNormalizedSH_sympyGeneratedExpressions(e,a)).sum()




if __name__ == "__main__" :
	sys.exit(unittest.main())

