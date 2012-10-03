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
from PySide.QtCore import Qt
import sys
import numpy as np
import time
from OpenGL import GL, GLU
import math

from colorfield import ColorField, Reloader

from audio3d.sphericalHarmonics import ead2xyz, semiNormalizedSH, sh, shi_reverse, on3d as orthonormalization, shGrid
from blobview import BlobView

def cartesian_product(*arrays):
	import operator
	broadcastable = np.ix_(*arrays)
	broadcasted = np.broadcast_arrays(*broadcastable)
	rows, cols = reduce(operator.mul, broadcasted[0].shape), len(broadcasted)
	out = np.empty(rows * cols, dtype=broadcasted[0].dtype)
	start, end = 0, rows
	for a in broadcasted:
		out[start:end] = a.reshape(-1)
		start, end = end, end + rows
	return out.reshape(cols, rows).T

def imageData(filename, width=None, height=None) :
	qimage = QtGui.QImage(filename)
	scaled = qimage if width is None and height is None else (
		qimage.scaled(width,height) )
	buffer = scaled.bits()
	return np.ndarray(
		shape = (height,width),
		dtype = np.uint8,
		buffer = buffer,
		)[:,:width].copy() # copy needed because the buffer is not persistent

def projectImageToSH(image) :
	"""
	Maps an array representing samples of a function defined in the surface
	of an sphere into spherical harmonics of the specified representation
	and order.
	Points are sampled at equal angular intervals for azimuth and elevation.
	(plate-carree projection).
	"""
	h,w = image.shape
	elevations, azimuths, sh = shGrid(h,w)
	cosines = np.cos(np.radians(elevations)).reshape(-1,1,1,1)
	sh = sh * cosines
	return image.reshape(w*h).dot(sh.reshape(w*h,-1)).reshape(sh.shape[2:])


class SphericalHarmonicsKnobs(QtGui.QWidget) :

	def __init__(self) :
		QtGui.QWidget.__init__(self)
		self._editing = True
		self._grid = QtGui.QGridLayout()
		self.setLayout(self._grid)

		def componentKnob(i,j) :
			l,m = shi_reverse(i,j)
			knob = QtGui.QDial()
			knob.setMinimum(-1000)
			knob.setMaximum(+1000)
			knob.setStyleSheet("background-color: %s"%orderColors[l%len(orderColors)])
			label = QtGui.QLabel("%d,%+d"%(l,m))
			label.setAlignment(Qt.AlignCenter)
			self._grid.addWidget(label, i*2, j)
			self._grid.addWidget(knob, i*2+1, j)
			knob.valueChanged.connect(self.knobEdited)
			return knob

		order = 5
		orderColors = [
			"#889988",
			"#667766",
			"#bb7788",
			"#aa6666",
			]
		self._knobs = [[
			componentKnob(i,j)
			for j in xrange(order+1) ]
			for i in xrange(order+1) ]


		self._editing = False

	functionChanged = QtCore.Signal()

	def sphericalHarmonicsMatrix(self) :
		return np.array([[
			knob.value()/1000.
			for knob in row ]
			for row in self._knobs ])

	def setSphericalHarmonicsMatrix(self, array) :
		self._editing = True
		array *= 1000./max(1.,abs(array).max())
		for i, row in enumerate(self._knobs) :
			for j,knob in enumerate(row) :
				knob.setValue(array[i,j])
		self._editing = False
		self.functionChanged.emit()

	def knobEdited(self) :
		if self._editing : return
		self.functionChanged.emit()

	def reset(self) :
		self._editing = True
		for row in self._knobs :
			for knob in row:
				knob.setValue(0)
		self._editing = False
		self.functionChanged.emit()

	def negate(self) :
		self._editing = True
		for row in self._knobs :
			for knob in row:
				knob.setValue(-knob.value())
		self._editing = False
		self.functionChanged.emit()



class SphereLab(QtGui.QWidget) :

	def __init__(self) :

		def addSpin(name, minimum, default, maximum, slot) :
			topLayout.addWidget(QtGui.QLabel(name+":"))
			spin = QtGui.QSpinBox()
			spin.setMinimum(minimum)
			spin.setMaximum(maximum)
			spin.setValue(default)
			spin.valueChanged.connect(slot)
			topLayout.addWidget(spin)
			return spin

		def addButton(layout, name, slot) :
			button = QtGui.QPushButton(name)
			button.clicked.connect(slot)
			layout.addWidget(button)


		QtGui.QWidget.__init__(self)
		self.shProjections = np.array([[[[]]]])
		self._editing = False
		self.setLayout(QtGui.QHBoxLayout())
		leftPanel = QtGui.QVBoxLayout()
		self.layout().addLayout(leftPanel)
		self._shknobs = SphericalHarmonicsKnobs()
		topLayout = QtGui.QHBoxLayout()
		self._parallelsSpin = addSpin("Parallels", 4, 50, 400, self.updateResolution)
		topLayout.addStretch(1)
		self._meridiansSpin = addSpin("Meridians", 4, 80, 400, self.updateResolution)

		addButton(topLayout, "Reset", self._shknobs.reset)
		addButton(topLayout, "Negate", self._shknobs.negate)
		addButton(topLayout, "Resynth", self.resynthesize)
		presetLayout1 = QtGui.QHBoxLayout()
		presetLayout2 = QtGui.QHBoxLayout()
		addButton(presetLayout1, "Front", self.sample_frontPoint)
		addButton(presetLayout1, "Back", self.sample_backPoint)
		addButton(presetLayout1, "Top", self.sample_topPoint)
		addButton(presetLayout1, "Down", self.sample_downPoint)
		addButton(presetLayout1, "Left", self.sample_leftPoint)
		addButton(presetLayout1, "Right", self.sample_rightPoint)
		addButton(presetLayout2, "Omni", self.sample_omni)
		addButton(presetLayout2, "Equator", self.sample_equator)
		addButton(presetLayout2, "Greenwitch", self.sample_greenwitch)
		addButton(presetLayout2, "Map", self.sample_map)

		leftPanel.addLayout(topLayout)
		leftPanel.addLayout(presetLayout1)
		leftPanel.addLayout(presetLayout2)
		leftPanel.addWidget(self._shknobs)


		rightPanel = QtGui.QVBoxLayout()
		self.layout().addLayout(rightPanel)

		self.blobView = BlobView()
		rightPanel.addWidget(self.blobView)

		self.synthetizedFunction = ColorField(width, height)
		rightPanel.addWidget(self.synthetizedFunction)

		self.targetFunction = ColorField(width, height)
		rightPanel.addWidget(self.targetFunction)

		rightPanel.setStretch(0,3)
		rightPanel.setStretch(1,1)
		rightPanel.setStretch(2,1)
		self.layout().setStretch(0,1)
		self.layout().setStretch(1,1)
		self._shknobs.functionChanged.connect(self.reloadData)

	def updateResolution(self) :
		self.reloadData()

	def sphericalHarmonicsMatrix(self) :
		return self._shknobs.sphericalHarmonicsMatrix()

	def setSphericalHarmonicsMatrix(self, array) :
		self._shknobs.setSphericalHarmonicsMatrix(array)

	def loadFromSamples(self, image) :
		h,w = image.shape
		print "Display data..."
		self.targetFunction.format(w, h, ColorField.signedScale)
		self.targetFunction.data()[:] = 127 + image*127./(max(1.,abs(image.max())))
		self.targetFunction.reload()

		print "Project to SH..."
		imageInSH = projectImageToSH(image)
		imageInSH *= orthonormalization

		print "Upload SH components to knobs..."
		self.setSphericalHarmonicsMatrix(imageInSH)
		print "Reloading..."
		self.reloadData()

	def reloadData(self) :
		nelevations = self._parallelsSpin.value()
		nazimuths = self._meridiansSpin.value()
		if self.shProjections.shape[:2] != (nelevations, nazimuths) :
			print "Reshaping %ix%x..."%(nelevations, nazimuths)
			self.elevations, self.azimuths, self.shProjections = shGrid(nelevations, nazimuths)
			self.shProjections *= orthonormalization
			print "Reshape outputs..."
			self.sphericalPoints = np.array([
				[self.elevations[ei], self.azimuths[ai], 0 ]
				for ei in xrange(nelevations)
				for ai in xrange(nazimuths)
				])
			self.indexes = np.array(
				[[
					[i+nazimuths*j,i+nazimuths*(j+1)]
					for i in xrange(nazimuths) ]
					for j in xrange(nelevations-1) ]
				).flatten()

		# taking coeficients from the knobs
		shMatrix = self.sphericalHarmonicsMatrix()

		self.data = self.shProjections.reshape(nazimuths*nelevations, shMatrix.size).dot(
			shMatrix.reshape(shMatrix.size )
			).reshape(self.shProjections.shape[:2])

		self.sphericalPoints[:,2] = (5*self.data).reshape(nelevations*nazimuths)

		self.blobView.setEadPoints(self.sphericalPoints)
		xyzs = np.array([ead2xyz(e,a,abs(d)) for e,a,d in self.sphericalPoints])
		self.blobView.scene()._vertices = xyzs
		self.blobView.scene()._normals = xyzs
		self.blobView.scene()._meshColors = np.array([
			[1.,.0,.0, .6] if d<0 else [0.,0.,1., .9]
			for e,a,d in self.sphericalPoints])
		self.blobView.scene()._indexes = self.indexes
		self.blobView.update()

		maxValue = abs(self.data).max()
		if maxValue > 1 : self.data /= maxValue
		self.synthetizedFunction.format(nazimuths, nelevations, ColorField.signedScale)
		self.synthetizedFunction.data()[:] = self.data/(2/255.)+127
		self.synthetizedFunction.reload()


	# Sphere sampling resolution for the synthetic data sets
	sampleResolution = 36*4, 18*4

	def resynthesize(self) :
		w,h = self.sampleResolution
		image = self.data
		self.loadFromSamples(image)

	def sample_map(self) :
		w,h = self.sampleResolution
		image = imageData("16bit_world_height.png", w, h)
		max = float(image.max())
		min = float(image.min())
		image = (image-min)/(max-min)
		self.loadFromSamples(image)

	def sample_omni(self) :
		w,h = self.sampleResolution
		image = np.ones((h,w))
		self.loadFromSamples(image)

	def sample_frontPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,w/2] = +1
		self.loadFromSamples(image)

	def sample_frontNegativePoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,w/2] = -1
		self.loadFromSamples(image)

	def sample_backNegativePoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,0] = -1
		self.loadFromSamples(image)

	def sample_backPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,0] = +1
		self.loadFromSamples(image)

	def sample_rightPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,w/4] = +1
		self.loadFromSamples(image)

	def sample_leftPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,3*w/4] = +1
		self.loadFromSamples(image)

	def sample_topPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[1,:] = +1 # north pole
		self.loadFromSamples(image)

	def sample_downPoint(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[-2,:] = +1 # south pole
		self.loadFromSamples(image)

	def sample_equator(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[h/2,:] = +1 # equator
		self.loadFromSamples(image)

	def sample_greenwitch(self) :
		w,h = self.sampleResolution
		image = np.zeros((h,w))
		image[:,w/2] = +1 # Grenwitch
		self.loadFromSamples(image)



if __name__ == "__main__" :

	width = 800
	height = 600

	app = QtGui.QApplication(sys.argv)

	w = QtGui.QDialog()
	w.setLayout(QtGui.QHBoxLayout())

	w0 = SphereLab()
	w.layout().addWidget(w0)
#	w.resize(width, height)

	from audio3d.sphericalHarmonics import shSize, shShape

	imageInSH = np.arange(shSize).reshape(shShape)*1000./shSize

	w0.setSphericalHarmonicsMatrix(imageInSH)

	w0.reloadData()

	w.show()

	app.exec_()



