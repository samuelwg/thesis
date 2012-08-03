#!/usr/bin/python

# This script shows graphs on how the spatial sampling rate is related
# to the accuracy at higher frequencies of the real Hrtf(t,theta) function.
# It takes the Spherical Harmonics Decomposition integrated using
# different stride steps.

import Gnuplot
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
import bmaudio
import numpy

configurations = [36,72]
configurations = [2,4,6,8,12,18,24,36,72]

import sys
component = "x"
if '--component' in sys.argv :
	component = sys.argv[sys.argv.index('--component')+1]

useLenght=201
#   y = e ^ -(s^2 w^2 / 2)
#   ln x = - 1/2 s^2 w^2 
#   -2 ln x = -s^2 w^2
#   sqrt(2 ln x) = s w

class SpectralPlot :
	def __init__(self,title=None):
		self.data=[]
		self.gp=Gnuplot.Gnuplot(persist=1)
		if title : self.gp.title(title)
		self.gp('set data style lines')
		self.gp.xlabel("Frequency (Hz)")
		self.gp.ylabel("Magnitude (dB)")
	def addSpectrum(self, spectrum, spectralRange, name) :
		self.addDbMagnitude(20*numpy.log10(abs(spectrum)), spectralRange, name)
	def addDbMagnitude(self, dbMagnitude, spectralRange, name) :
		nBins = len(dbMagnitude)
		self.data.append( Gnuplot.Data (
			spectralRange/nBins*numpy.arange(nBins), 
			dbMagnitude,
			title=name
			))
	def show(self) :
#		self.gp("set yrange [-100:20]")
		self.gp("set xrange [0:24000]")
		self.gp.plot(*self.data)
	def hardcopy(self, file, format='png') :
		print "Storing '%s'"%file
		self.gp.hardcopy(filename=file,terminal=format) 

correlations = dict()
divergenceThreshold = .005
errors = dict()

def equivalentFile(database, speakers, component) :
	dboptions = dict(
		kreuzer="--kreuzer",
		mitFull="",
		mitDiffuse="--mitdiffuse",
		)
	databasefile = bmaudio.selectHrtfDatabase( [dboptions[database], "--speakers", str(speakers) ])
	return bmaudio.hrtfDatabaseToEquivalentPath(databasefile, component)

for database in [
	"kreuzer",
	"mitFull",
	"mitDiffuse",
	] :
	plot = SpectralPlot("Divergence from a 72 speakers configuration (%s component, %s database)"%(component,database))
	sr,spectrum0 = bmaudio.fileToSpectrum(equivalentFile(database,72,component))
	if database == "kreuzer" :
		spectrum0 = spectrum0[:201]
	errors[database] = []
	correlations[database] = []
	spectroPlotData = []
	for configuration in configurations :
		sr,spectrum = bmaudio.fileToSpectrum(equivalentFile(database,configuration,component))
		if database == "kreuzer" :
			spectrum = spectrum[:201]
			sr = 40200
		diffs = abs(spectrum-spectrum0)#/abs(spectrum0)
#		spectroPlotData.append(20*numpy.log10(diffs))
		minDiff = diffs[:]
		minDiff[minDiff<1e-10] = 1e-4
		spectroPlotData.append(20*numpy.log10(minDiff))
#		spectroPlotData.append(diffs)
		cumsum = diffs.cumsum()
		plot.addSpectrum(diffs,sr/2,"%s speakers"%configuration)
#		plot.addSpectrum(cumsum,sr/2,"%s speakers"%configuration)
		nonZeros = numpy.where(diffs>divergenceThreshold)[0]
		nonZeros = numpy.where(cumsum>.05)[0]
		correlations[database].append(nonZeros[0]*sr/2/len(spectrum0) if len(nonZeros) else sr/2)
		errors[database].append(numpy.sqrt(sum(abs(spectrum-spectrum0)**2))/len(spectrum))
	plot.show()
	plot.hardcopy("%s2D-divergence.eps"%database,"postscript")
	import pylab
	pylab.pcolor(numpy.array(spectroPlotData),cmap=pylab.cm.bone)
	pylab.yticks(range(len(configurations)), configurations)
	pylab.title(database)
	pylab.show()

print correlations
gp=Gnuplot.Gnuplot(persist=1)
gp.title("Frequencies at which divergence with %i configuration is greater than %f"%(configurations[-1],divergenceThreshold))
gp.ylabel("Frequency (Hz)")
gp.xlabel("Number of speakers")
gp('set data style linespoints')
gp.plot(*[
	Gnuplot.Data(
		configurations[:-1],
		correlation[:-1],
		title=database)
	for database, correlation in correlations.iteritems()])
gp.hardcopy(filename="correlation.eps", terminal="postscript")

print errors
gp=Gnuplot.Gnuplot(persist=1)
gp.title("Quadratic mean difference with the %i configuration"%(configurations[-1]))
gp.ylabel("Quadratic mean error")
gp.xlabel("Number of speakers")
gp('set data style linespoints')
#gp('set log y')
gp.plot(*[
	Gnuplot.Data(
		configurations[:-1],
		data[:-1],
		title=database)
	for database, data in errors.iteritems()])
gp.hardcopy(filename="correlation.eps", terminal="postscript")


