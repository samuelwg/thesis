#!/usr/bin/env python

"""
Plotting routines to encapsulate frequent plot configuration
"""
import numpy
import pylab
import bmaudio
import math


def plotField(basename, field, xaxis, yaxis,
		levels=numpy.arange(-60,6,3),
		minxdb=-60,
		logFreq = False,
		) :
	field[numpy.isinf(field)]=-6
	field-=field.max()
	pylab.rcParams['font.size'] = 14
	pylab.rcParams['figure.figsize']=(15,6)
	pylab.rcParams['figure.subplot.left'] = 0.06
	pylab.rcParams['figure.subplot.right'] = 0.99
	pylab.contourf(xaxis, yaxis, field, levels, cmap=pylab.cm.Greys_r, extend="both")
	pylab.colorbar(
		fraction=.06,
		pad=0.02,
		).set_label("Magnitude (dB)")
	c=pylab.contour(xaxis, yaxis, field, levels, cmap=pylab.cm.Greys)
	pylab.ylim(-90,90)
	pylab.yticks([-90,-45,0,+45,+90])
	if logFreq :
		pylab.semilogx()
		pylab.xticks(
			[1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000],
			[1,2,5,10,20,50,100,200,500,'1k','2k','5k','10k','20k'],
			)
		pylab.xlim(10) # start with 10^1
	else :
		pylab.xticks(
			[i*1000 for i in xrange (0,23,2)],
			["%ik"%i  if i else "0" for i in xrange (0,23,2)],
			)
		
	pylab.xlabel("Frequency (Hz)")
	pylab.ylabel("Source azimuth (degrees)")
	pylab.clabel(c, fmt="%i")
	pylab.savefig(basename+".pdf", format="pdf")
	pylab.savefig(basename+".png", format="png")
	pylab.show()
	pylab.close()


class SpectrumDisplay :
	def __init__(self) :
		self.data = []
		self.phase = []
		self._inDb = False
		self._showPhase = False
		self.ymin = None
		self.ymax = None
		self.fmin = None
		self.fmax = None
		self._logFreq = False
		self.stylePreference = ['colors','lines']
		self.styleVariation = dict(
			colors = "rgbcmyk",
			lines = ["-", "--", ":", "-."],
			)
		self.legendTitle = None
		self.colorTitle = None
		self.colorItems = None
		self.lineTitle = None
		self.lineItems = None
		self._legendPosition = "lower right"

	def legendPosition(self, position) :
		self._legendPosition = position
	def setLegendTitle(self, title) :
		self.legendTitle=title
	def setStyleLegend(self, colorTitle, colorItems, lineTitle, lineItems) :
		self.colorTitle = colorTitle
		self.colorItems = colorItems
		self.lineTitle = lineTitle
		self.lineItems = lineItems
	def inDb(self) :
		self._inDb = True
	def setStyleVariation(self, name, values) :
		self.styleVariation[name] = values[:]
	def setStylePreference(self, stylePreference) :
		self.stylePreference = stylePreference
	def ylim(self, ymin=None, ymax=None) :
		self.ymin = ymin
		self.ymax = ymax
	def flim(self, fmin=None, fmax=None) : 
		self.fmin = fmin
		self.fmax = fmax
	def showPhase(self, showPhase=True) :
		self._showPhase = showPhase;
	def setLogFrequency(self) :
		self._logFreq=True
	def addWaveFile(self, file) :
		samplerate, spectrum = fileToSpectrum(file)
		self.addSpectrumData(spectrum, samplerate/2, file)
	def addSpectrumData(self, spectrum, spectralRange, name, style=None) :
		nBins = len(spectrum)
		if self._inDb :
			magnitude = 20*numpy.log10(bmaudio.removeZeros(abs(spectrum)))
		else :
			magnitude = abs(spectrum)
		frequencies =  numpy.array(xrange(nBins))*spectralRange/nBins
		self.data.append( ( frequencies, magnitude, name))
		self.phase.append( ( frequencies, numpy.angle(spectrum), name))

	def _style(self, i, styleVars) :
		result= ""
		previousStyle = []
		for styleVar in styleVars :
			result += styleVar[i%len(styleVar)]
			i//=len(styleVar)
		return result

	def _build(self) :
		import pylab
		styleVars = [self.styleVariation[style] for style in self.stylePreference]
		pylab.rcParams['figure.figsize']=(15,6)
		pylab.rcParams['figure.subplot.left'] = 0.05
		pylab.rcParams['figure.subplot.right'] = 0.95
		pylab.rcParams['figure.subplot.bottom'] = 0.09
		pylab.rcParams['figure.subplot.top'] = 0.95
		pylab.clf()
		if self._showPhase :
			pylab.subplot(2,1,1)
		if self._logFreq :
			pylab.semilogx()
			pylab.xticks(
				[1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000],
				[1,2,5,10,20,50,100,200,500,'1k','2k','5k','10k','20k'],
				)
			pylab.xlim(10) # start with 10^1
		else :
			pylab.xticks(
				[i*1000 for i in xrange (0,23,2)],
				["%ik"%i  if i else "0" for i in xrange (0,23,2)],
				)
		if self._inDb :
			pylab.ylabel("Magnitude (dB)")
		else:
			pylab.ylabel("Magnitude")
		pylab.xlabel("Frequency (Hz)")
		useLabels = self.colorTitle or self.lineTitle
		for i, (frequencies, values, name) in enumerate(self.data) :
			pylab.plot(frequencies, values, self._style(i, styleVars), label=None if useLabels else name)
		colorLegend = None
		lineLegend = None
		if self.colorTitle :
			colorArtists = [
				pylab.Rectangle((0,0),1,1, fc=color)
				for color in self.styleVariation['colors'][:len(self.colorItems)]
				]
			colorLegend = pylab.legend(colorArtists, self.colorItems, loc="upper right", title=self.colorTitle)
			lineArtists = [
				pylab.Line2D([1],[1], linestyle=lineStyle, color="b")
				for lineStyle in self.styleVariation['lines'][:len(self.lineItems)]
				]
			lineLegend = pylab.legend(lineArtists, self.lineItems, loc="lower right", title=self.lineTitle)
		elif self.legendTitle :
			pylab.legend(loc=self._legendPosition, title=self.legendTitle)
		else:
			pylab.legend(loc=self._legendPosition)
		pylab.axis([self.fmin, self.fmax, self.ymin, self.ymax])
		ax = pylab.gca()
		print "fmax: ", self.fmax
		ax.set_autoscale_on(False)

		pylab.ylim(self.ymin, self.ymax)
		pylab.xlim(self.fmin, self.fmax)
		if self._showPhase :
			pylab.subplot(2,1,2)
			for i, (frequencies, values, name) in enumerate(self.phase) :
				pylab.plot(frequencies, values, self._style(i, styleVars), label=None if useLabels else name)
			pylab.ylabel("Phase (Radians)")
			pylab.xlabel("Frequency")
			pylab.ylim(-math.pi, math.pi)
			if self._logFreq :
				pylab.semilogx()
				pylab.xticks(
					[1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000],
					[1,2,5,10,20,50,100,200,500,'1k','2k','5k','10k','20k'],
					)
				if not self.fmin : self.fmin = 10 # start with 10^1 if it was 0
			else :
				pylab.xticks(
					[i*1000 for i in xrange (0,23,2)],
					["%ik"%i  if i else "0" for i in xrange (0,23,2)],
					)
			pylab.xlim(self.fmin, self.fmax)
		if colorLegend : pylab.gca().add_artist(colorLegend)
		if lineLegend : pylab.gca().add_artist(lineLegend)

	def show(self) :
		import pylab
		self._build()
		pylab.show()
	def hardcopy(self, file, format='png') :
		import pylab
		self._build()
		pylab.savefig(file,format=format, papersize='a4')
		pylab.show()

def displayArray( audioArray, samplerate, name, inDb=True ):
	displaySpectrum(numpy.fft.rfft(audioArray), samplerate/2, name, inDb)

def displaySpectrum( spectrum, spectralRange, name, inDb=True, showPhase=False ):
	Sd = SpectrumDisplay()
	if inDb: Sd.inDb()
	if showPhase: Sd.showPhase()
	Sd.addSpectrumData(spectrum, spectralRange, name)
	Sd.show()

def displayWaveFile( file, name, inDb=True ):
	Sd = SpectrumDisplay()
	Sd.addWaveFile( file )
	if inDb:
		Sd.inDb()
	Sd.show()




