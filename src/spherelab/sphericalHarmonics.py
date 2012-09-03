#!/usr/bin/python

from PySide import QtCore, QtGui, QtOpenGL
import sys
import numpy as np
import time
import math


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

def sh(sh, e, a) :
	x,y,z = ead2xyz(e, a, 1)
	return math.sqrt(1/math.pi)*(
		math.sqrt(1./2) * (
			sh[0,0] +
			0
		) +
		math.sqrt(1./4) * (
			sh[0,1] * x +
			sh[1,1] * z +
			sh[1,0] * y +
			0
		) +
		(
			sh[2,0] * (x*y)             * math.sqrt(15./4) +
			sh[2,1] * (x*z)             * math.sqrt(15./4) +
			sh[2,2] * (2*z*z -x*x -y*y) * math.sqrt( 5./16) +
			sh[1,2] * (y*z)             * math.sqrt(15./4) +
			sh[0,2] * (x*x -y*y)        * math.sqrt(15./16) +
			0
		) +
		(
			sh[3,0] * x*(x*x-3*y*y)           * math.sqrt(35./32) +
			sh[3,1] * z*x*y                   * math.sqrt(105./4)+
			sh[3,2] * x*(4*z*z -x*x -y*y)     * math.sqrt(21./32) +
			sh[3,3] * z*(2*z*z -3*x*x -3*y*y) * math.sqrt(7./16) +
			sh[2,3] * y*(4*z*z -x*x -y*y)     * math.sqrt(21./32) +
			sh[1,3] * z*(x*x-y*y)             * math.sqrt(105./16) +
			sh[0,3] * y*(3*x*x-y*y)           * math.sqrt(35./32) +
			0
		) +
		0
	)



def shIndexes(order) :
	return [(l,m) for l in xrange(order+1) for m in xrange(-l,l+1) ]

def shIndex2Matrix(l,m) :
	return (l-m,l) if m>0 else (l,l+m)

import unittest

class SphericalHarmonicsTests(unittest.TestCase) :
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

		self.assertEqual((0,0), shIndex2Matrix(0,0))

		self.assertEqual((1,0), shIndex2Matrix(1,-1))
		self.assertEqual((1,1), shIndex2Matrix(1,0))
		self.assertEqual((0,1), shIndex2Matrix(1,1))
                             
		self.assertEqual((2,0), shIndex2Matrix(2,-2))
		self.assertEqual((2,1), shIndex2Matrix(2,-1))
		self.assertEqual((2,2), shIndex2Matrix(2,0))
		self.assertEqual((1,2), shIndex2Matrix(2,1))
		self.assertEqual((0,2), shIndex2Matrix(2,2))


	def test_sh_0_0(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(0,0)] = 1
		self.assertAlmostEqual(1, sh(components, 90, 0)*math.sqrt(2*math.pi))
		self.assertAlmostEqual(1, sh(components, -90, 0)*math.sqrt(2*math.pi))
		self.assertAlmostEqual(1, sh(components, 30, 4)*math.sqrt(2*math.pi))

	def test_sh_1_1(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(1,+1)] = 1
		max = math.sqrt(1./4/math.pi)
		self.assertAlmostEqual(sh(components,   0,   0)/max,+1)
		self.assertAlmostEqual(sh(components,   0, 180)/max,-1)
		self.assertAlmostEqual(sh(components, -90,   0)/max, 0)
		self.assertAlmostEqual(sh(components, +90,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, -90)/max, 0)
		self.assertAlmostEqual(sh(components,   0, +90)/max, 0)
		self.assertAlmostEqual(sh(components,  30,   4)/max, 0.86391580942710433)

	def test_sh_1_m1(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(1,-1)] = 1
		max = math.sqrt(1./4/math.pi)
		self.assertAlmostEqual(sh(components,   0,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, 180)/max, 0)
		self.assertAlmostEqual(sh(components, -90,   0)/max, 0)
		self.assertAlmostEqual(sh(components, +90,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, -90)/max,-1)
		self.assertAlmostEqual(sh(components,   0, +90)/max,+1)
		self.assertAlmostEqual(sh(components,  30,   4)/max, 0.060410878340834702)

	def test_sh_1_0(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(1,0)] = 1
		max = math.sqrt(1./4/math.pi)
		self.assertAlmostEqual(sh(components,   0,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, 180)/max, 0)
		self.assertAlmostEqual(sh(components, -90,   0)/max,-1)
		self.assertAlmostEqual(sh(components, +90,   0)/max,+1)
		self.assertAlmostEqual(sh(components,   0, -90)/max, 0)
		self.assertAlmostEqual(sh(components,   0, +90)/max, 0)
		self.assertAlmostEqual(sh(components,  30,   0)/max, 0.5)

	def test_sh_2_2(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(2,+2)] = 1
		max = math.sqrt(15./16/math.pi)
		self.assertAlmostEqual(sh(components,   0,   0)/max,+1)
		self.assertAlmostEqual(sh(components,   0, 180)/max,+1)
		self.assertAlmostEqual(sh(components, -90,   0)/max, 0)
		self.assertAlmostEqual(sh(components, +90,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, -90)/max,-1)
		self.assertAlmostEqual(sh(components,   0, +90)/max,-1)
		self.assertAlmostEqual(sh(components,  30,   4)/max, 0.7427010515561776)

	def test_sh_2_1(self) :
		components = np.zeros((4,4))
		components[shIndex2Matrix(2,+1)] = 1
		max = math.sqrt(15./16/math.pi)
		self.assertAlmostEqual(sh(components,   0,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, 180)/max, 0)
		self.assertAlmostEqual(sh(components, -90,   0)/max, 0)
		self.assertAlmostEqual(sh(components, +90,   0)/max, 0)
		self.assertAlmostEqual(sh(components,   0, -90)/max, 0)
		self.assertAlmostEqual(sh(components,   0, +90)/max, 0)
		self.assertAlmostEqual(sh(components,   0,  45)/max, 1)
		self.assertAlmostEqual(sh(components,   0, -45)/max, 1)
		self.assertAlmostEqual(sh(components,   0, 135)/max, 1)
		self.assertAlmostEqual(sh(components,   0, 315)/max, 1)
		self.assertAlmostEqual(sh(components,   0, -90)/max, 0)
		self.assertAlmostEqual(sh(components,   0, +90)/max, 0)
		self.assertAlmostEqual(sh(components,  30,   4)/max, 0.7427010515561776)


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

