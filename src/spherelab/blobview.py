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
from OpenGL import GL, GLU
import math

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


class BlobViewScene(QtGui.QGraphicsScene) :

	def __init__(self) :
		super(BlobViewScene, self).__init__()
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

		super(BlobViewScene, self).drawBackground(painter,rect)
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



class BlobView(QtGui.QGraphicsView) :

	def __init__(self) :
		super(BlobView, self).__init__()
		self.setRenderHints(QtGui.QPainter.Antialiasing | QtGui.QPainter.SmoothPixmapTransform)
#		self.setRenderHints(QtGui.QPainter.SmoothPixmapTransform)
		viewPort = QtOpenGL.QGLWidget()
		self.setViewport(viewPort)
		self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate);
		scene = BlobViewScene()
		self.setScene(scene);


	def resizeEvent(self, event) :
		if self.scene() :
			self.scene().setSceneRect(QtCore.QRect(QtCore.QPoint(0, 0), event.size()))
		super(BlobView, self).resizeEvent(event)

	def setEadPoints(self, points) :
		self.scene().setEadPoints( points)



if __name__ == "__main__" :

	width = 800
	height = 600
	import sys
	import numpy as np

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



