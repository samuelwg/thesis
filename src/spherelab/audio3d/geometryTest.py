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


if __name__ == "__main__":
	unittest.main()


