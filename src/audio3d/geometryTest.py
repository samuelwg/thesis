#!/usr/bin/python

import unittest

from geometry import *

class AngleNormalizationTest(unittest.TestCase):

	def setUp(self) :
		pass

	def tearDown(self):
		pass

	def testNormalizeAngle_withZeroAngles(self) :
		self.assertEquals(
			(0.0, 0.0),
			normalizeAngles(0.0, 0.0))

	def testNormalizeAngle_inPositiveQuadrant(self) :
		self.assertEquals(
			(10.0, 20.0),
			normalizeAngles(10.0, 20.0))

	def testNormalizeAngle_withNegativeAzimut(self) :
		self.assertEquals(
			(350.0, 10.0),
			normalizeAngles(-10.0, 10.0))

	def testNormalizeAngle_withAzimutBeyond360(self) :
		self.assertEquals(
			(10.0, 10.0),
			normalizeAngles(370.0, 10.0))

	def testNormalizeAngle_withElevationBeyond90(self) :
		self.assertEquals(
			(190.0, 89.0),
			normalizeAngles(10.0, 91.0))


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


if __name__ == "__main__":
	unittest.main()


