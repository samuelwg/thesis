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
			[round(v,6) for v in expected],
			[round(v,6) for v in result],
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

	def test_xyz2ead_front(self) :
		self.assertCoordsEqual( xyz2ead(*[1.0, 0.0, 0.0]), (0, 0, 1) )

	def test_xyz2ead_further_front(self) :
		self.assertCoordsEqual( xyz2ead(*[2.0, 0.0, 0.0]), (0, 0, 2) )

	def test_xyz2ead_back_turning_left(self) :
		self.assertCoordsEqual( xyz2ead(*[-1.0, 0.0, 0.0]), (0, 180, 1) )

	def test_xyz2ead_left(self) :
		self.assertCoordsEqual( xyz2ead(*[0.0, 1.0, 0.0]), (0, 90, 1) )

	def test_xyz2ead_right(self) :
		self.assertCoordsEqual( xyz2ead(*[0.0, -1.0, 0.0]), (0, -90, 1) )

	def test_xyz2ead_up(self) :
		self.assertCoordsEqual( xyz2ead(*[0.0, 0.0, 1.0]), (90, 0, 1) )

	def test_xyz2ead_down(self) :
		self.assertCoordsEqual( xyz2ead(*[0.0, 0.0, -1.0]), (-90, 0, 1) )

	def test_xyz2ead_mid_elevated(self) :
		self.assertCoordsEqual( xyz2ead(*[0.8660254, 0.0, 0.5]), (30, 0, 1) )

	def test_xyz2ead_mixed(self) :
		self.assertCoordsEqual( xyz2ead(*[0.75, 0.4330127, 0.5]), (30, 30, 1) )

class ChordDistanceTests(unittest.TestCase) :

	def test_chordDistance_equal(self) :
		self.assertEqual(0,  chordDistance(30, 3, 30, 3) )

	def test_chordDistance_antipodal(self) :
		self.assertEqual(4,  chordDistance(30, 0, -30, +180) )

	def test_chordDistance_quarter(self) :
		self.assertAlmostEqual(2, chordDistance(0, 0, 0, +90) )

	def test_chordDistance_movingAlongEquatorLeft(self) :
		self.assertAlmostEqual(0.26794919243112258,  chordDistance(0, 0, 0, +30) )

	def test_chordDistance_movingAlongEquatorRight(self) :
		self.assertAlmostEqual(0.26794919243112258,  chordDistance(0, 0, 0, -30) )

	def test_chordDistance_movingAlongGreenwichUp(self) :
		self.assertAlmostEqual(0.26794919243112258,  chordDistance(0, 0, +30, 0) )

	def test_chordDistance_movingAlongGreenwichDown(self) :
		self.assertAlmostEqual(0.26794919243112258,  chordDistance(0, 0, -30, 0) )

	def test_chordDistance_movingAlongGreenwichUpper(self) :
		self.assertAlmostEqual(0.26794919243112258,  chordDistance(30, 0, 60, 0) )

	def test_chordDistance_azimuthAngleAreSmallerDistanceInUpperElevation(self) :
		self.assertAlmostEqual(0.20096189432334199,  chordDistance(30, 0, 30, +30) )








if __name__ == "__main__":
	unittest.main()


