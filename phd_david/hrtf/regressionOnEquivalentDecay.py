#!/usr/bin/python
import Gnuplot
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
import bmaudio
import numpy
import math

configurations = [36,72]
configurations = [
	2,4,6,8,
	12,18,24,36,72]
component = "w"
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
		self.gp("set yrange [-80:30]")
		self.gp.plot(*self.data)
	def hardcopy(self, file, format='png') :
		print "Storing '%s'"%file
		self.gp.hardcopy(filename=file,terminal=format) 

correlations = dict()
divergenceThreshold = .005


slopes = {}
pearsonRs = {}

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
	
	baseName = "interleaved/%s2D-%%02i-E%%s.wav"%database
	plot = SpectralPlot("(%s component, %s database)"%(component,database))
	databaseSlope = []
	pearsonR = []
	for configuration in configurations :
		sr,spectrum = bmaudio.fileToSpectrum(equivalentFile(database,configuration,component))
		if database == "kreuzer" :
			spectrum = spectrum[:200]
			sr = 40000
		freqs = sr/2./len(spectrum)*numpy.arange(len(spectrum))
		x = (freqs**2)
		y = numpy.log(abs(spectrum))
		slope = sum(x * y)/sum(x**2)
		regression = numpy.exp(slope*(freqs**2))
#		plot.addSpectrum(y,sr/2,"%s speakers"%configuration)
		plot.addSpectrum(regression,sr/2,"%s speakers"%configuration)
		c=343
		xx = sum(x**2)
		xy = sum(x*y)
		yy = sum(y**2)
		y2 = sum(y)**2
		n = len(y)
		ym = sum(y)/n # y.mean()
		m = xy/xx
#		r2 = 1 - sum((y-m*x)**2) / sum((y-ym))**2)
#		r2 = 1 - (yy -2*m*xy + m**2 * xx) / (yy -2*ym*sum(y) + n*ym*ym)
#		r2 = 1 - (yy -2*xy*xy/xx + xy*xy/xx) / (yy -2*ym*n*ym + n*ym*ym)
#		r2 = 1 - (yy -xy*xy/xx) / (yy -n*ym*ym)
#		r2 = (yy -n*yn*yn - yy +xy*xy/xx) / (yy -n*ym*ym)
#		r2 = (-n*yn*yn + xy*xy/xx) / (yy -n*ym*ym)
		r2 = (xy*xy/xx - y2/n ) / (yy - y2/n)
		databaseSlope.append(100*2*c*math.sqrt(-slope if slope<0 else 0)/math.pi)
#		pearsonR.append(sum( ((x-x.mean())/x.std()) * ((y-y.mean())/y.std()) ) / len(x))
		pearsonR.append(r2)
	slopes[database] = databaseSlope
	pearsonRs[database] = pearsonR
	plot.show()
	plot.hardcopy("%s2D-regression.eps"%database,"postscript")
print pearsonR
gp=Gnuplot.Gnuplot(persist=1)
gp.title("Slopes")
gp.ylabel("Standard deviation (cm)")
gp.xlabel("Number of speakers")
gp('set data style linespoints')
gp.plot(*[
	Gnuplot.Data(
		configurations,
		databaseSlope,
		title=database)
	for database, databaseSlope in slopes.iteritems()])
gp.hardcopy(filename="slopes.eps", terminal="postscript")


gp=Gnuplot.Gnuplot(persist=1)
gp.title("How well spectrum fits decay")
gp.ylabel("Pearson r")
gp.xlabel("Number of speakers")
gp('set data style linespoints')
gp.plot(*[
	Gnuplot.Data(
		configurations,
		r,
		title=database)
	for database, r in pearsonRs.iteritems()])
gp.hardcopy(filename="slopes.eps", terminal="postscript")




