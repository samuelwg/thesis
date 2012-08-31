#!/usr/bin/python

from PySide import QtCore, QtGui
import numpy as np
import time


class ColorField(QtGui.QWidget) :

	greyScale = [QtGui.qRgb(i,i,i) for i in xrange(256)]
	fancyScale = [
		QtGui.qRgb(i,0,0) for i in xrange(0,256,4) ] + [
		QtGui.qRgb(0xff-i,i,0) for i in xrange(0,256,4) ] + [
		QtGui.qRgb(0,0xff-i,i) for i in xrange(0,256,4) ] + [
		QtGui.qRgb(i,i,0xff) for i in xrange(0,256,4) ]
	signedScale = [
		QtGui.qRgb(255-2*i,12,12) for i in xrange(0,125) ] + [
		QtGui.qRgb(2,0x2f,0)  ] + [
		QtGui.qRgb(1,0x6f,0)  ] + [
		QtGui.qRgb(0,0xff,0)  ] + [
		QtGui.qRgb(0,0x6f,1)  ] + [
		QtGui.qRgb(0,0x2f,2)  ] + [
		QtGui.qRgb(12,12,2*(i-130)) for i in xrange(130,256) ]

	def __init__(self, width, height, colorTable = None) :
		super(ColorField, self).__init__()
		self.resize(width, height)
		self.format(width, height, colorTable)

	def format(self, width, height, colorTable = None) :
		self._colorTable = colorTable
		qformat, npformat, pixelBytes = (
				QtGui.QImage.Format_Indexed8, np.uint8, 1
			) if colorTable else (
				QtGui.QImage.Format_RGB32, np.uint32, 4
			)

		self._img = QtGui.QImage(width, height, qformat)
		if colorTable : self._img.setColorTable(colorTable)

		self._rawdata = np.ndarray(
			shape=(height,self._img.bytesPerLine()/pixelBytes),
			dtype=npformat,
			buffer=self._img.bits())[:,:width]


	def data(self) :
		return self._rawdata

	def paintEvent(self, event) :
		painter = QtGui.QPainter(self)
		painter.drawImage(self.rect(), self._img)

	def reload(self) :
		self.update()


class Reloader(QtCore.QObject) :
	def __init__(self, target) :
		super(Reloader,self).__init__()
		self._target = target
		self._counter = 0
		self._tic = time.clock()

	def timerEvent(self, event) :
		data = self._target.data()
		data[:] = np.random.randint(0,0x00ffffff,size=data.shape).astype(data.dtype) # 12m
#		data[:] = 0x00ff0000
		data[self._counter:self._counter+1,10:-10] = 0x0000ff00
#		data[:256,-4] = np.arange(256)
#		data[:256,-3] = np.arange(256)
#		data[:256,-2] = np.arange(256)
#		data[:256,-1] = np.arange(256)
		self._counter += 1
		self._target.reload()
		period = 200
		if not self._counter%period :
			toc = time.clock()
			print "%.2f frames/s"% (period/float(toc - self._tic) if toc!=self._tic else 0), "%2.2f ms"% (float(toc - self._tic)*1000/period)
			self._tic = toc

if __name__ == "__main__" :

	import sys
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

