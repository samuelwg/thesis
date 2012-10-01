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

from audio3d.sphericalHarmonics import ead2xyz, semiNormalizedSH, sh, shi_reverse

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
	cosines = np.cos(np.radians(elevations))
	sh = sh*cosines.reshape(-1,1,1,1)
	shShape = sh.shape[2:]
	shSize = shShape[0]*shShape[1]
	# normalize, regarding its concentration on higher elevations
#	for ei, e in enumerate(elevations) :
#		sh[ei] = np.cos(np.radians(e))

	result = image.reshape(w*h).dot(sh.reshape(w*h,(shSize))).reshape(shShape)
	return result


class memoize:
	def __init__(self, function):
		self.function = function
		self.memoized = {}
	def __call__(self, *args):
		try:
			return self.memoized[args]
		except KeyError:
			self.memoized[args] = self.function(*args)
			return self.memoized[args]


@memoize
def shGrid(nelevations, nazimuths) :
	"""Given a number of elevation and azimuth sampling position returns
	the elevations and azimuths and a matrix ne X na X sho X sho
	the sampling of the Seminormalized Spherical Harmonic functions
	for the cross product of the elevation and the azimuth sampling.
	"""
	# inverted elevation order
	elevations = np.linspace(90,  -90, nelevations)
	azimuths = np.linspace( -180, 180, nazimuths, endpoint=False)
	shsampling = np.array([[
		semiNormalizedSH(e,a)
		for a in azimuths]
		for e in elevations]
		)
	return elevations, azimuths, shsampling


class TrackBall(object) :
	def __init__(self, angularVelocity=0., axis=QtGui.QVector3D(0,1,0)) :
		self._angularVelocity = angularVelocity
		self._rotation = QtGui.QQuaternion()
		self._paused = False
		self._pressed = False
		self._lastTime = QtCore.QTime.currentTime()
		self._axis = axis

	def rotation(self) :
		if self._paused or self._pressed : return self._rotation
		currentTime = QtCore.QTime.currentTime()
		angle = self._angularVelocity * self._lastTime.msecsTo(currentTime)
		return QtGui.QQuaternion.fromAxisAndAngle(self._axis, angle) * self._rotation

	def push(self, point, reference) :
		self._rotation = self.rotation()
		self._pressed = True
		self._lastTime = QtCore.QTime.currentTime()
		self._lastPos = point
		self._angularVelocity = 0

	def move(self, point, reference) :
		if not self._pressed : return

		currentTime = QtCore.QTime.currentTime()
		msecs = self._lastTime.msecsTo(currentTime)
		if msecs <= 20 : return # ignore frequenta

		def map3d(point2) :
			point3 = QtGui.QVector3D(point2)
			sqrZ = 1 - point2.x()**2 - point2.y()**2
			if sqrZ > 0 :
				point3.setZ(math.sqrt(sqrZ))
			else :
				point3.normalize()
			return point3

		dotProduct = QtGui.QVector3D.dotProduct
		crossProduct = QtGui.QVector3D.crossProduct

		if False : # Plane method
			delta = QtCore.QLineF(self._lastPos, point)
			self._angularVelocity = math.degrees(delta.length()/msecs)
			self._axis = QtGui.QVector3D(-delta.dy(), delta.dx(), 0.0).normalized()
			self._axis = reference.rotatedVector(self._axis)
			self._rotation = QtGui.QQuaternion.fromAxisAndAngle(self._axis, math.degrees(delta.length())) * self._rotation
		else : # Sphere method

			lastPos3D = map3d(self._lastPos)
			currentPos3D = map3d(point)

			self._axis = crossProduct(lastPos3D, currentPos3D)
			angle = math.degrees(math.asin(math.sqrt(dotProduct(self._axis, self._axis))))

			self._angularVelocity = angle / msecs
			self._axis.normalize()
			self._axis = reference.rotatedVector(self._axis)
			self._rotation = QtGui.QQuaternion.fromAxisAndAngle(self._axis, angle) * self._rotation

		self._lastTime = currentTime
		self._lastPos = point


	def release(self, point, reference) :
		self.move(point, reference)
		self._pressed = False


class SpherePointScene(QtGui.QGraphicsScene) :

	def __init__(self) :
		super(SpherePointScene, self).__init__()
		self._frame = 0

		self._trackballs = [
			TrackBall(),
			TrackBall(),
			TrackBall(),
			]
		self.startTimer(20)
		self._points = []
		self._vertices = None
		self._indexes = None
		self._normals = None
		self._distance = 20

	def timerEvent(self, event) :
		self.update()

	def setEadPoints(self, points) :
		self._points = points
		self.update()

	def drawBackground(self, painter, rect) :
		height = float(painter.device().height())
		width = float(painter.device().width())

		painter.beginNativePainting()
		self.setStates()

		GL.glClearColor(0.0, 0.0, 0.0, 0.0)
		GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

		GL.glMatrixMode(GL.GL_PROJECTION)
		GLU.gluPerspective(60.0, width / height, 0.01, 450.0);
		GLU.gluLookAt(-self._distance,0,0,0,0,0,0,0,1);

		GL.glMatrixMode(GL.GL_MODELVIEW);

		view = QtGui.QMatrix4x4()
		view.rotate(self._trackballs[2].rotation())

#		view = np.array(list(view.data())).reshape((4,4))
#		view[2, 3] -= 2.0 * math.exp(self._distance / 1200.0)
#		view = QtGui.QMatrix4x4(*view.reshape((16,)))

		GL.glLoadMatrixf(view.data())
		self.drawAxis()
#		self.drawPoints()
		self.drawSphere()

		self.setDefaultState()
		painter.endNativePainting()
		self._frame+=1

	def drawPoints(self) :
		GL.glEnable(GL.GL_CULL_FACE)
		GL.glEnable(GL.GL_LIGHTING)
		quadric = GLU.gluNewQuadric()
		for e,a,d in self._points :
			ra, re = math.radians(a), math.radians(e)
			sa, se = math.sin(ra), math.sin(re)
			ca, ce = math.cos(ra), math.cos(re)
			if d<0 :
				GL.glColor(1., 1., 1., .5)
				d = -d
			else :
				GL.glColor(1., .6, .6, .5)
			GL.glPushMatrix()
			GL.glTranslate(d*ce*ca, d*ce*sa,d*se)
			GLU.gluSphere(quadric, .1, 10, 10)
			GL.glPopMatrix()

	def drawSphere(self) :
		if self._indexes is None : return

		GL.glDisable(GL.GL_CULL_FACE)
#		GL.glEnable(GL.GL_LIGHTING)
#		GL.glScale(1,1,1)

#		GL.glDisable(GL.GL_LIGHTING)
		GL.glEnable(GL.GL_BLEND)
		GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)

		GL.glEnableClientState( GL.GL_VERTEX_ARRAY )
		GL.glEnableClientState( GL.GL_COLOR_ARRAY )
		GL.glEnableClientState( GL.GL_NORMAL_ARRAY )

		GL.glColorPointer( 4, GL.GL_FLOAT, 0, self._meshColors)
		GL.glVertexPointer( 3, GL.GL_FLOAT, 0, self._vertices )
		GL.glNormalPointer( GL.GL_FLOAT, 0, self._normals )
		GL.glDrawElements( GL.GL_TRIANGLE_STRIP, len(self._indexes), GL.GL_UNSIGNED_INT, self._indexes )

		GL.glDisableClientState( GL.GL_COLOR_ARRAY )
		GL.glColor(0.6,.3,.6)
		GL.glDrawElements( GL.GL_POINTS, len(self._indexes), GL.GL_UNSIGNED_INT, self._indexes )

		GL.glDisableClientState( GL.GL_COLOR_ARRAY )
		GL.glDisableClientState( GL.GL_VERTEX_ARRAY )
		GL.glDisableClientState( GL.GL_NORMAL_ARRAY )


	def drawAxis(self) :

		GL.glPushAttrib(GL.GL_ENABLE_BIT)
		GL.glDisable(GL.GL_LIGHTING)

		GL.glLineWidth(1)

		GL.glBegin(GL.GL_TRIANGLE_FAN)
		GL.glColor(0,0,1.)
		GL.glVertex3f(0, 0, +10)
		GL.glVertex3f(-.3, +.3, +9.5)
		GL.glVertex3f(+.3, +.3, +9.5)
		GL.glVertex3f(+.3, -.3, +9.5)
		GL.glVertex3f(-.3, -.3, +9.5)
		GL.glVertex3f(-.3, +.3, +9.5)
		GL.glEnd()

		GL.glBegin(GL.GL_TRIANGLE_FAN)
		GL.glColor(0,1.,0)
		GL.glVertex3f(0, +10, 0)
		GL.glVertex3f(-.3, +9.5, +.3)
		GL.glVertex3f(+.3, +9.5, +.3)
		GL.glVertex3f(+.3, +9.5, -.3)
		GL.glVertex3f(-.3, +9.5, -.3)
		GL.glVertex3f(-.3, +9.5, +.3)
		GL.glEnd()

		GL.glBegin(GL.GL_TRIANGLE_FAN)
		GL.glColor(1.,0,0)
		GL.glVertex3f(+10, 0, 0)
		GL.glVertex3f(+9.5, -.3, +.3)
		GL.glVertex3f(+9.5, +.3, +.3)
		GL.glVertex3f(+9.5, +.3, -.3)
		GL.glVertex3f(+9.5, -.3, -.3)
		GL.glVertex3f(+9.5, -.3, +.3)
		GL.glEnd()

		GL.glBegin(GL.GL_LINES)

		GL.glColor(0,0,1.)
		GL.glVertex3f(0, 0, -10)
		GL.glVertex3f(0, 0, +10)

		GL.glColor(1.,0,0)
		GL.glVertex3f(-10, 0, 0)
		GL.glVertex3f(+10, 0, 0)

		GL.glColor(0,1.,0)
		GL.glVertex3f(0, -10, 0)
		GL.glVertex3f(0, +10, 0)

		GL.glEnd()

		GL.glPopAttrib()


	def setStates(self) :
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_CULL_FACE)
		GL.glEnable(GL.GL_LIGHTING)
		GL.glEnable(GL.GL_COLOR_MATERIAL)
#		GL.glEnable(GL.GL_TEXTURE_2D)
		GL.glEnable(GL.GL_NORMALIZE)
		GL.glPolygonMode(GL.GL_FRONT, GL.GL_LINE);

		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glPushMatrix()
		GL.glLoadIdentity()

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glPushMatrix()
		GL.glLoadIdentity()

		self.setLights()

		materialSpecular = [0.2, 0.5, 0.5, 1.0]
		GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, materialSpecular);
		GL.glMaterialf(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, 32),

	def setDefaultState(self) :
#		GL.glClearColor(0.0, 0.0, 0.0, 0.0)

		GL.glDisable(GL.GL_DEPTH_TEST)
		GL.glDisable(GL.GL_CULL_FACE)
		GL.glDisable(GL.GL_LIGHTING)
#		GL.glDisable(GL.GL_COLOR_MATERIAL)
		GL.glDisable(GL.GL_TEXTURE_2D)
		GL.glDisable(GL.GL_LIGHT0)
		GL.glDisable(GL.GL_NORMALIZE)

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glPopMatrix()

		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glPopMatrix()

		GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, 0.0)
		defaultMaterialSpecular = [0.0, 0.0, 0.0, 1.0]
		GL.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, defaultMaterialSpecular)
		GL.glMaterialf(GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, 0.0)

	def setLights(self) :

		GL.glColorMaterial(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT_AND_DIFFUSE)
		lightColour = 1.0, 1.0, 0.8, 0.1
		lightDiffuse = .4, .4, 0.4, 1.
		lightDir = 0.0, 1.0, 0.8, 0.0
		GL.glLightfv(GL.GL_LIGHT0, GL.GL_DIFFUSE, lightDiffuse)
		GL.glLightfv(GL.GL_LIGHT0, GL.GL_SPECULAR, lightColour)
		GL.glLightfv(GL.GL_LIGHT0, GL.GL_POSITION, lightDir);
		GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0)
		GL.glEnable(GL.GL_LIGHT0)

		GL.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE, (.4,.5,.7))
		GL.glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, lightColour)
		GL.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, (-.5,-.2,-.5,0));
		GL.glLightModelf(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0)
		GL.glEnable(GL.GL_LIGHT1)


	def wheelEvent(self, event) :
		QtGui.QGraphicsScene.wheelEvent(self,event)
		if event.isAccepted() : return

		self._distance += event.delta()/120.
		if self._distance < 10  : self._distance = 10
		if self._distance > 50  : self._distance = 50
		event.accept()

	def pixelPosToViewPos(self, p) :
		return QtCore.QPointF(
			1.0 - 2.0 * float(p.y()) / self.height(),
			2.0 * float(p.x()) / self.width() - 1.0,
			)

	def mouseMoveEvent(self, event) :
		QtGui.QGraphicsScene.mouseMoveEvent(self, event)
		if event.isAccepted() : return

		mousePos = self.pixelPosToViewPos(event.scenePos())

		if event.buttons() & Qt.LeftButton :
			self._trackballs[0].move(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()
		else :
			self._trackballs[0].release(mousePos,
				self._trackballs[2].rotation().conjugate())

		if event.buttons() & Qt.RightButton :
			self._trackballs[1].move(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()
		else :
			self._trackballs[1].release(mousePos,
				self._trackballs[2].rotation().conjugate());


		if event.buttons() & Qt.MidButton :
			self._trackballs[2].move(mousePos, QtGui.QQuaternion())
			event.accept();
		else :
			self._trackballs[2].release(mousePos, QtGui.QQuaternion())



	def mousePressEvent(self, event) :
		QtGui.QGraphicsScene.mousePressEvent(self, event)
		if event.isAccepted() : return

		mousePos = self.pixelPosToViewPos(event.scenePos())

		if event.buttons() & Qt.LeftButton :
			self._trackballs[0].push(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()

		if event.buttons() & Qt.RightButton :
			self._trackballs[1].push(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()

		if event.buttons() & Qt.MidButton :
			self._trackballs[2].push(mousePos, QtGui.QQuaternion())
			event.accept();


	def mouseReleaseEvent(self, event) :
		QtGui.QGraphicsScene.mouseReleaseEvent(self, event)
		if event.isAccepted() : return

		mousePos = self.pixelPosToViewPos(event.scenePos())

		if event.buttons() & Qt.LeftButton :
			self._trackballs[0].release(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()

		if event.buttons() & Qt.RightButton :
			self._trackballs[1].release(mousePos,
				self._trackballs[2].rotation().conjugate())
			event.accept()

		if event.buttons() & Qt.MidButton :
			self._trackballs[2].release(mousePos, QtGui.QQuaternion())
			event.accept();



class SpherePointView(QtGui.QGraphicsView) :

	def __init__(self) :
		super(SpherePointView, self).__init__()
		self.setRenderHints(QtGui.QPainter.Antialiasing | QtGui.QPainter.SmoothPixmapTransform)
#		self.setRenderHints(QtGui.QPainter.SmoothPixmapTransform)
		viewPort = QtOpenGL.QGLWidget()
		self.setViewport(viewPort)
		self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate);
		scene = SpherePointScene()
		self.setScene(scene);


	def resizeEvent(self, event) :
		if self.scene() :
			self.scene().setSceneRect(QtCore.QRect(QtCore.QPoint(0, 0), event.size()))
		super(SpherePointView, self).resizeEvent(event)

	def setEadPoints(self, points) :
		self.scene().setEadPoints( points)


class SphereLab(QtGui.QWidget) :

	def __init__(self) :
		self.shProjections = np.array([[[[]]]])

		def addSpin(name, minimum, default, maximum, slot) :
			topLayout.addWidget(QtGui.QLabel(name+":"))
			spin = QtGui.QSpinBox()
			spin.setMinimum(minimum)
			spin.setMaximum(maximum)
			spin.setValue(default)
			spin.valueChanged.connect(slot)
			topLayout.addWidget(spin)
			return spin

		QtGui.QWidget.__init__(self)
		self._editing = False
		self.setLayout(QtGui.QHBoxLayout())
		leftPanel = QtGui.QVBoxLayout()
		self.layout().addLayout(leftPanel)
		self._grid = QtGui.QGridLayout()
		topLayout = QtGui.QHBoxLayout()
		self._parallelsSpin = addSpin("Parallels", 4, 50, 400, self.resolutionChanged)
		topLayout.addStretch(1)
		self._meridiansSpin = addSpin("Meridians", 4, 80, 400, self.resolutionChanged)

		def addButton(layout, name, slot) :
			button = QtGui.QPushButton(name)
			button.clicked.connect(slot)
			layout.addWidget(button)

		addButton(topLayout, "Reset", self.reset)
		addButton(topLayout, "Negate", self.negate)
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
		leftPanel.addLayout(self._grid)


		def componentKnob(i,j) :
			l,m = shi_reverse(i,j)
			knob = QtGui.QDial()
			knob.setMinimum(-100)
			knob.setMaximum(+100)
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

		rightPanel = QtGui.QVBoxLayout()
		self.layout().addLayout(rightPanel)

		self.blobView = SpherePointView()
		rightPanel.addWidget(self.blobView)

		self.synthetizedFunction = ColorField(width, height)
	#	reloader1 = Reloader(self.synthetizedFunction)
	#	reloader1.startTimer(0)
		rightPanel.addWidget(self.synthetizedFunction)

		self.targetFunction = ColorField(width, height)
		rightPanel.addWidget(self.targetFunction)

		rightPanel.setStretch(0,3)
		rightPanel.setStretch(1,1)
		rightPanel.setStretch(2,1)
		self.layout().setStretch(0,1)
		self.layout().setStretch(1,1)
		self.resolutionChanged.connect(self.reloadData)
		self.functionChanged.connect(self.reloadData)



	functionChanged = QtCore.Signal()
	resolutionChanged = QtCore.Signal()

	def knobEdited(self) :
		if self._editing : return
		self.functionChanged.emit()

	def sphericalHarmonicsMatrix(self) :
		return np.array([[
			knob.value()/100.
			for knob in row ]
			for row in self._knobs ])

	def setSphericalHarmonicsMatrix(self, array) :
		self._editing = True
		for i, row in enumerate(self._knobs) :
			for j,knob in enumerate(row) :
				knob.setValue(array[i,j])
		self._editing = False
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

	def loadFromSamples(self, image) :
		h,w = image.shape
		imageBytesIn8Bits = 127 + image*127./(max(1.,abs(image.max())))

		self.targetFunction.format(w, h, ColorField.signedScale)
		self.targetFunction.data()[:,:] = imageBytesIn8Bits
		self.targetFunction.reload()

		imageInSH = projectImageToSH(image)
		imageInSH *= 100./max(1.,abs(imageInSH).max())

		self.setSphericalHarmonicsMatrix(imageInSH)
		self.reloadData()

	def reloadData(self) :
		nelevations = self._parallelsSpin.value()
		nazimuths = self._meridiansSpin.value()
		if self.shProjections.shape[:2] != (nelevations, nazimuths) :
			print "Reshaping..."
			self.elevations, self.azimuths, self.shProjections = shGrid(nelevations, nazimuths)
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

		data = self.shProjections.reshape(nazimuths*nelevations, shMatrix.size).dot(
			shMatrix.reshape(shMatrix.size )
			).reshape(self.shProjections.shape[:2])

		self.sphericalPoints[:,2] = (5*data).reshape(nelevations*nazimuths)

		self.blobView.setEadPoints(self.sphericalPoints)
		xyzs = np.array([ead2xyz(e,a,abs(d)) for e,a,d in self.sphericalPoints])
		self.blobView.scene()._vertices = xyzs
		self.blobView.scene()._normals = xyzs
		self.blobView.scene()._meshColors = np.array([[1.,.0,.0, .6] if d<0 else [0.,0.,1., .9] for e,a,d in self.sphericalPoints])
		self.blobView.scene()._indexes = self.indexes
		self.blobView.update()
		maxValue = abs(data).max()
		if maxValue > 1 : data /= maxValue
		self.synthetizedFunction.format(nazimuths, nelevations, ColorField.signedScale)
		self.synthetizedFunction.data()[:] = data/(2/255.)+127
		self.synthetizedFunction.reload()


	# Sphere sampling resolution for the synthetic data sets
	sampleResolution = 36*4, 18*4

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

	imageInSH = np.arange(shSize).reshape(shShape)*100./shSize

	w0.setSphericalHarmonicsMatrix(imageInSH)

	w0.reloadData()

	w.show()

	app.exec_()



