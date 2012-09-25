#!/usr/bin/env python

"""
This module contains common tools for the audio team of Barcelona Media.
"""

import wave, numpy, os, math
import sys
import decorator
from audio3d.plot import SpectrumDisplay

@decorator.decorator
def verbose(f, self, *args, **keyargs) :
	"This method decorator allows to trace a call to the decorated function"
	print "[%i] %s.%s( %s )" % (
		id(self),
		self.__class__.__name__,
		f.__name__,
		", ".join(
			[ str(arg) for arg in args ] +
			[name+"="+str(arg) for name,arg in keyargs.iteritems()])
		)
	result = f(self, *args, **keyargs)
	if result is not None : print "=>", result
	return result


# do not use this
_hrtf_path = None
def hrtf_path() :
	# Lazy initialization for _hrtf_path.
	# it avoids forcing having HRTF_PATH defined if it is not needed.
	global _hrtf_path
	if _hrtf_path : return _hrtf_path
	if 'HRTF_PATH' not in os.environ.keys():
		print >> sys.stderr, 'Could not find HRTF_PATH environment variable. Please define it and point it to where imm_bm/HRTFs is'
		sys.exit(-1)
	_hrtf_path = os.environ['HRTF_PATH']
	return _hrtf_path


def saveWave(filename, data, samplerate, verbose=False) :
	if verbose: print "Writting %s..."%filename
	wavefile = wave.open( filename, 'w' )
	wavefile.setframerate( samplerate )
	wavefile.setnchannels( 1 )
	wavefile.setsampwidth( 2 ) # number of bytes per sample
	wavefile.writeframesraw((32767.0*data).astype(numpy.int16).tostring())
	wavefile.close()

def loadWave(filename, verbose=False) :
	if verbose: print "Loading %s..."%filename
	wavefile = wave.open( filename, 'r' )
	samplingRate   = wavefile.getframerate()
	bytesPerSample = wavefile.getsampwidth()
	nFrames        = wavefile.getnframes()
	data           = wavefile.readframes( nFrames )
	wavefile.close()
	if bytesPerSample==2:  # 16 bits per sample
		data = numpy.fromstring( data, numpy.int16 ) / 32767.0 # -1..1 values
	elif bytesPerSample==1: # 8 bits per sample
		data = (numpy.fromstring( data, numpy.uint8 ) / 128.0 ) - 1.0 # -1..1 values
	elif bytesPerSample==4 : # 32 bits per sample
		data = (numpy.fromstring( data, numpy.int32 ) / 2.**(32-1) ) # -1..1 values
	else :
		print bytesPerSample
		2+"asdfasdf"
	return samplingRate, data

def firstMaxima(data, threshold=.85) :
	reference = threshold * max(abs(data))
	for i, (v0, v1, v2) in  enumerate(zip(data[:-2],data[1:-1],data[2:])) :
		if abs(v1) <= reference : continue
		if v0>v1 and v1<v2: return i+interpolateMaxima(v0,v1,v2)[0] # maximum
		if v0<v1 and v1>v2: return i+interpolateMaxima(v0,v1,v2)[0] # minimum
	return 0

def interpolateMaxima(y0,y1,y2) :
	"""Given three equidistant values for y so that y1 is either greater 
	or less than both y0 and y2, this function returns the interpolated
	maxima or minima (xm,ym). It considers that x0=0, x1=1 and x2=2.
	Multiply xm by the real x1'-x0' and add x0' to get the real x position.
	The formula uses the derivative of the quadratic lagrange interpolation.
	The formula rationale is in CLAM::CircularPeakPeaker::interpolate.
	"""
	a = y0/2 - y1 + y2/2
	b = y1 -y0 -a
	c = y0
	xmax = -b/(a*2)
	ymax = b*xmax/2 + y0
	return xmax, ymax

# warning: this code is duplicated in scripts/AngleNormalizationTest.py,
# update it as well!!!
def normalizeAngles(azimuth, elevation) :
	elevation += 90.
	elevation %= 360.
	if elevation > 180. :
		azimuth += 180.
		elevation = 360 - elevation
	elevation -= 90.
	azimuth %= 360.
	return azimuth,elevation

# algorithm copied from livecoreo/widget/Geometry.hxx
def spherical2Cartesian(azimuthDegrees, elevationDegrees, distance) :
	a = math.radians(azimuthDegrees)
	e = math.radians(elevationDegrees)
	ce = math.cos(e)
	ca = math.cos(a)
	se = math.sin(e)
	sa = math.sin(a)
	return distance * ce * ca, distance * ce * sa, distance * se

def cartesian2Spherical(x, y, z) :
	"""Returns spherical coordinates for x,y,z in a tuple azimuth, elevation, distance"""

	d = math.sqrt(x*x+y*y+z*z)

	if (d < 0.00001):
		return 0, 0, 0

	return math.degrees(math.atan2(y,x)), math.degrees(math.asin(z/d)), d

def angularDistance(a1,e1,a2,e2):
	""" Implements the haversine formula (see Wikipedia) """
	return 2*math.asin(math.sqrt((math.sin((e1 - e2)/2))**2 + math.cos(e1)*math.cos(e2)*(math.sin((a1 - a2)/2))**2));

#algorithm copied from AbsoluteCoordinates2RelativeAngles.hxx (CLAM/plugins/spacialization)
def absoluteCoordinates2RelativeAngles(listenerX,listenerY,listenerZ,listenerAzimuthDegrees,listenerElevationDegrees,listenerRollDegrees,sourceX,sourceY,sourceZ):
	listenerAzimuth=math.radians(listenerAzimuthDegrees)
	listenerElevation=math.radians(listenerElevationDegrees)
	listenerRoll=math.radians(listenerRollDegrees)		
	dx=(sourceX-listenerX)
	dy=(sourceY-listenerY)
	dz=(sourceZ-listenerZ)
	cosAzimuth=math.cos(listenerAzimuth)
	sinAzimuth=math.sin(listenerAzimuth)
	cosElevation=math.cos(listenerElevation)
	sinElevation=math.sin(listenerElevation)
	cosRoll=math.cos(listenerRoll)
	sinRoll=math.sin(listenerRoll)
	cosASinE=cosAzimuth*sinElevation
	sinASinE=sinAzimuth*sinElevation
	rotatedX = dx*cosAzimuth*cosElevation + dy*sinAzimuth*cosElevation + dz*sinElevation
	rotatedY = dx*(-sinAzimuth*cosRoll -cosASinE*sinRoll) + dy*(cosAzimuth*cosRoll - sinASinE*sinRoll) + dz*(cosElevation*sinRoll)
	rotatedZ = dx*(-cosASinE*cosRoll + sinAzimuth*sinRoll) + dy*(-sinASinE*cosRoll - cosAzimuth*sinRoll) + dz*(cosElevation*cosRoll)
	azimuth=math.atan2(rotatedY,rotatedX)
	module=math.sqrt(rotatedX*rotatedX+rotatedY*rotatedY+rotatedZ*rotatedZ)
	if module < 1e-12:
		elevation=rotatedZ
	else:
		elevation=math.asin(rotatedZ/module)
	return math.degrees(azimuth), math.degrees(elevation), module

def calculateGainByDistance(distance,exponent=1.0,minimumDistance=1.0,distanceThreshold=0.0):
	gain=0
	if minimumDistance==0:
		minimumDistance=0.0001  #avoid zero div.
	if distance<minimumDistance:
		distance=minimumDistance
	if distanceThreshold==0 or distance<=distanceThreshold:
		gain=1.0/math.pow(distance,exponent)
	return gain

def dbToAmplitude(db):
	return pow(10, (db/20.))

def dbToEnergy(db):
	return pow(10, (db/10.))

def amplitudeToDb(amplitude):
	return 20 * math.log(amplitude, 10.)

def energyToDb(energy):
	return 10 * math.log(energy, 10.)

class HrtfDatabase :

	def __init__(self, databaseFile) :
		self._databaseFile = databaseFile
		self._data = []
		self._audio = {}
		base = os.path.dirname(databaseFile)
		for line in open(databaseFile) :
			try : elevation, azimuth, filename = line.split()
			except: continue
			azimuth, elevation = normalizeAngles(float(azimuth), float(elevation))
			self._data.append( ( float(elevation), float(azimuth), os.path.join(base,filename)) )
		self._orientationToFilename = dict(((e,a),f) for e,a,f in self._data)
		self._data.sort()

	def equivalentPath(self,component) :
		return hrtfDatabaseToEquivalentPath(self._databaseFile, component)

	def layout(self) :
		layout = {}
		for elevation, azimuth, filename in self._data :
			layout[elevation] = layout.get(elevation, 0) + 1
		return layout

	def __repr__(self) :
		return "%s:%r" % ( self.__class__.__name__, self._data )

	def printOrientations(self) :
		print "printOrientations", self.layout()
		for elevation, n in self.layout() :
			print elevation, []

	def lowerElevation(self) :
		"""Returns the lower elevation"""
		return min( elevation for elevation, azimuth, filename in self._data )

	def nearestOrientation(self, e, a):
		return min((angularDistance(a, e, azimuth, elevation),(elevation, azimuth)) for elevation, azimuth, filename in self._data)[1]

	def preload(self):
		for elevation, azimuth, filename in self._data:
			loadHrtf(elevation, azimuth)

	def loadHrtf(self, elevation, azimuth) :
		samplerate, wave = loadWave(self._orientationToFilename[elevation,azimuth])
		self._audio[elevation,azimuth] = wave
		return wave

	def hrtfWave(self, elevation, azimuth) :
		if (elevation,azimuth) not in self._orientationToFilename :
			elevation, azimuth = self.nearestOrientation(elevation, azimuth)
		if (elevation, azimuth) in self._audio : return self._audio[elevation,azimuth]
		return self.loadHrtf(elevation, azimuth)


	def capCompensation(self, lowerElevation) :
		"""
		Given an speaker layout with horizontal symmetry but missing speakers on the bottom cap,
		returns the number of missing speakers and two factors for the pressure and the
		velocity of the lower elevation speakers in order to compensate the missing speakers.
		"""
		layout = self.layout()
		if len(layout) == 1 : return 0, 1., 1.
		lowerElevationSize = layout[lowerElevation]
		print lowerElevationSize
		# Includes the missing and the edge
		missingAngles = [ (float(-angle),float(n)) for angle,n in layout.iteritems() if -angle<=lowerElevation ]
		vCompensation = sum( nSpeakers * math.sin(math.radians(angle)) for angle,nSpeakers in missingAngles)
		vCompensation/= lowerElevationSize * math.sin(math.radians(lowerElevation))
		pCompensation = sum(nSpeakers for angle,nSpeakers in missingAngles) / lowerElevationSize
		nMissingAngles = sum( n for angle, n in layout.iteritems() if -angle not in layout.keys())
		print "Lower elevation:", lowerElevation, "nMissingAngles:", nMissingAngles, "pCompensation:", pCompensation, "vCompensation:", vCompensation
		return nMissingAngles, pCompensation, vCompensation

def elevationAndAzimuthFromName(filename, pattern) :
	"""\
	Extracts the elevation and azimuth values from the name.
	It uses a regular expression with named matching groups 
	'elevation' and 'azimuth'.
	It expects degrees on the name and returns radians.
	"""
	values = pattern.search(filename).groupdict()
	elevation = float(values['elevation'])*math.pi/180.
	while elevation > math.pi: elevation -= 2*math.pi
	azimuth = float(values['azimuth'])*math.pi/180.
	return elevation, azimuth

def fileToSpectrum(file) :
	samplerate, wave = loadWave(file)
	return samplerate, numpy.fft.rfft(wave)

def spectrumMagnitudeInDb(spectrum) :
	return 20*numpy.log(abs(spectrum))

def spectrumPhase(spectrum) :
	return numpy.angle(spectrum)

def minimumPhase(buffer) :
	zeroPaddedSize = len(buffer)*8
	return numpy.fft.irfft(minimumPhaseSpectrum(abs(numpy.fft.rfft(buffer,zeroPaddedSize))))[:len(buffer)]

def removeZeros(absoluteArray) :
	return numpy.array([1e-14 if val<1e-14 else val for val in absoluteArray])

def minimumPhaseSpectrum(spectrumMagnitude) :
	realCepstrum = numpy.fft.irfft(numpy.log(removeZeros(spectrumMagnitude)))
	halfLen = (len(realCepstrum)-1)/2
	oddLen  = (len(realCepstrum)-1)%2
	filter = numpy.array([1] + [2]*halfLen + [1]*oddLen + [0]*halfLen)
	complexSpectrum = filter*realCepstrum
	return numpy.exp(numpy.fft.rfft(complexSpectrum))

def delayWithMinimumPhaseCorrelation(buffer, filterFreq=None, interpolate=True, sampleRate=44100) :
	"""Returns the time of arrival of an impulse response.
	The computation is done by cross-correlating, in the frequency domain,
	the impulse response with its minimum phase equivalent.
	If filterFreq is given, those frequencies are removed from the original
	impulse response before all the computation.
	If interpolate is True then the delay is given with more precission than a sample
	by doing quadratic interpolation.
	"""
	spectrum = numpy.fft.rfft(buffer)
	if not filterFreq is None :
		spectrum[filterFreq*len(spectrum)/(sampleRate/2.):]=0
	undelayed = minimumPhaseSpectrum(numpy.abs(spectrum))
	correlation = spectrum*undelayed.conjugate()
	timeCorrelation = numpy.fft.irfft(correlation)
	maxDelay = numpy.argmax(timeCorrelation)
	if not interpolate : return maxDelay
	interpolatedX, interpolatedY = interpolateMaxima(
		timeCorrelation[maxDelay-1 if maxDelay>0 else len(timeCorrelation)-1],
		timeCorrelation[maxDelay],
		timeCorrelation[maxDelay+1 if maxDelay+1<len(timeCorrelation) else 0]
		)
	return maxDelay+interpolatedX-1

def delay85(buffer) :
	maximum = abs(buffer).max()
	return min([len(buffer)]+[i for i,x in enumerate(buffer) if abs(x)>maximum*.85])



def selectHrtfDatabase( args ) :

	if '--database' in args :
		return args[args.index('--database')+1]

	earLetter = 'R' if '-r' in args else 'L'
	speakerInfix = ''
	if '--speakers' in args :
		# TODO: Interleaved should work again
		density = int(args[args.index('--speakers')+1]) # TODO: Remove it
		speakerInfix="-2d%02i"%density
	if '--ircam' in args :
		# IRCAM database place left at 90 degrees
		nIrcam = int(args[args.index('--ircam')+1])
		return os.path.join(hrtf_path(),"ircam%i%s.hrtfs"%(nIrcam, earLetter))
	elif '--kreuzer' in args :
		# KREUZER database place left at 90 degrees
		return os.path.join(hrtf_path(),"kreuzerDatabase2%s.hrtfs"%(speakerInfix))
	else : # MIT KEMAR database
		# MIT database place left at 270 degrees
		if '--mitdiffuse' in args :
			return os.path.join(hrtf_path(),"mitKemarDiffuse%s%s.hrtfs" % (earLetter, speakerInfix))
		elif '--mitcompact' in args :
			return os.path.join(hrtf_path(),"mitKemarCompact%s%s.hrtfs" % (earLetter, speakerInfix))
		else :
			return os.path.join(hrtf_path(),"mitKemarFull%s%s.hrtfs" % (earLetter, speakerInfix))

def hrtfDatabaseToEquivalentPath(databaseFile, component) :
	return os.path.join(
		os.path.dirname(__file__), 
		"../../../bformat2binaural/equivalentIRs",
		os.path.splitext(os.path.split(databaseFile)[1])[0],
		"E%s.wav"%component,
		)

def headlessDisplacedDelaySimplified(azimuthRadians, headRadius, distance) :

	"""Returns the distance walked by a sound comming from azimuthRadians
	at given distance of the center of an spheric head with the given headRadius but ignoring any head.
	"""
	r = distance
	R = headRadius
	a = azimuthRadians-math.pi/2 # move the 0 angle to the right ear
	a0 = math.acos(R/r)
	"""
	while a<0 : a+=2*math.pi # Periodic each 2pi
	while a>2*math.pi : a-=2*math.pi # Periodic each 2pi
	if a>math.pi: a = 2*math.pi-a # Symmetric on pi on an spheric head
	"""
	return r - R*math.cos(a)



def sphericalHeadDelaySimplified(azimuthRadians, headRadius, distance) :
	"""Returns the distance walked by a sound comming from azimuthRadians
	at given distance of the center of an spheric head with the given headRadius.
	"""
	r = distance
	R = headRadius
	a = azimuthRadians-math.pi/2 # move the 0 angle to the right ear
	a0 = math.acos(R/r)
	while a<0 : a+=2*math.pi # Periodic each 2pi
	while a>2*math.pi : a-=2*math.pi # Periodic each 2pi
	if a>math.pi: a = 2*math.pi-a # Symmetric on pi on an spheric head
	return (
		r - R*math.cos(a)
			if a<math.pi/2 else
		r + R*(a - math.pi/2)
		)

def sphericalHeadDelay(azimuthRadians, headRadius, distance) :
	"""Returns the distance walked by a sound comming from azimuthRadians
	at given distance of the center of an spheric head with the given headRadius.
	Simplified by considering distance >> headRadius
	"""
	r = distance
	R = headRadius
	a = azimuthRadians-math.pi/2 # move the 0 angle to the right ear
	a0 = math.acos(R/r)
	while a<0 : a+=2*math.pi # Periodic each 2pi
	while a>2*math.pi : a-=2*math.pi # Periodic each 2pi
	if a>math.pi: a = 2*math.pi-a # Symmetric on pi on an spheric head
	return (
		math.sqrt(R*R + r*r -2*R*r*math.cos(a))
			if a<a0 else
		R*(a - a0) + math.sqrt(r*r - R*R)
		)

def azimuthDelaysFromHorizontalHrtf(hrtfDatabase) :
	"""Returns the samplingrate of the hrtfs and a dictionary of azimuthDegrees -> delay,
	where the delay is the time of arrival relative to the first arrival orientation."""
	delays = {}
	for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
		from math import radians
		if elevationDegrees!= 0 : continue
		samplingRate, data = loadWave(response)
		t = delayWithMinimumPhaseCorrelation(data, filterFreq=20000, interpolate=True) / samplingRate
		delays[azimuthDegrees] = t
	# Remove the minimal delay
	minDelay = min(delays.values())
	for azimuth, delay in delays.iteritems(): delays[azimuth] -= minDelay

	return delays

class SympyEvaluator :
	def __init__(self, function, *args) :
		self._variables = args
		self._function = function
	def __call__(self, *args) :
		return self._function.subs(dict(zip(self._variables,args))).evalf()
	def __repr__(self) :
		return "f(" + (",".join( [var.name for var in self._variables])) + ") = " + str(self._function)

def sphericalHarmonic(l,n) :
	"""Returns an evaluator for the seminormalized real spherical harmonic
	of order l, degree n as defined in:
	http://ambisonics.iem.at/xchange/format/ambisonics-xchange-format-appendix
	"""
	import sympy
	sign = -1 if n<0 else +1
	n = abs(n)
	x,y=sympy.var("x y")
	f= sympy.powsimp(sympy.trigsimp(
		sympy.assoc_legendre(l,n,sympy.sin(x)) *
		(-1)**n *
		sympy.sqrt((2*l+1) *
		(1 if n==0 else 2) *
		sympy.factorial(l-n)/sympy.factorial(l+n)) *
		(sympy.cos(n*y) if sign>=0 else sympy.sin(n*y))
		))
	return SympyEvaluator(f,x,y)

sphericalHarmonics2D = [ #  order, name, (m,n,rho), function, weight
	(0, 'w', (0,0,+1), lambda azimuth, elevation : 1, 1./numpy.sqrt(2)),
	(1, 'x', (1,1,+1), lambda azimuth, elevation : numpy.cos(azimuth), 1.),
	(1, 'y', (1,1,-1), lambda azimuth, elevation : numpy.sin(azimuth), 1.),
	(2, 'u', (2,2,+1), lambda azimuth, elevation : numpy.cos(2*azimuth), 2./numpy.sqrt(3)),
	(2, 'v', (2,2,-1), lambda azimuth, elevation : numpy.sin(2*azimuth), 2./numpy.sqrt(3)),
	(3, 'p', (3,3,+1), lambda azimuth, elevation : numpy.cos(3*azimuth), numpy.sqrt(8./5.)),
	(3, 'q', (3,3,-1), lambda azimuth, elevation : numpy.sin(3*azimuth), numpy.sqrt(8./5.)),
] # taken from http://www.york.ac.uk/inst/mustech/3d_audio/higher_order_ambisonics.pdf

sphericalHarmonics = [ #  order, name, (m,n,rho), function, weight
	(0, 'w', (0,0,+1), lambda azimuth, elevation : 1, 1./numpy.sqrt(2)),
	(1, 'z', (1,0,+1), lambda azimuth, elevation : numpy.sin(elevation), 1.),
	(1, 'x', (1,1,+1), lambda azimuth, elevation : numpy.cos(elevation)*numpy.cos(azimuth), 1.),
	(1, 'y', (1,1,-1), lambda azimuth, elevation : numpy.cos(elevation)*numpy.sin(azimuth), 1.),
	(2, 'r', (2,0,+1), lambda azimuth, elevation : (3*numpy.sin(elevation)**2-1)/2, 1.),
	(2, 's', (2,1,+1), lambda azimuth, elevation : numpy.sqrt(3./4)*numpy.cos(2*elevation)*numpy.cos(azimuth), 2./numpy.sqrt(3)),
	(2, 't', (2,1,-1), lambda azimuth, elevation : numpy.sqrt(3./4)*numpy.cos(2*elevation)*numpy.sin(azimuth), 2./numpy.sqrt(3)),
	(2, 'u', (2,2,+1), lambda azimuth, elevation : numpy.sqrt(3./4)*(numpy.cos(elevation)**2)*numpy.cos(2*azimuth), 2./numpy.sqrt(3)),
	(2, 'v', (2,2,-1), lambda azimuth, elevation : numpy.sqrt(3./4)*(numpy.cos(elevation)**2)*numpy.sin(2*azimuth), 2./numpy.sqrt(3)),
	(3, 'k', (3,0,+1), lambda azimuth, elevation : numpy.sin(elevation)*(5*numpy.sin(elevation)-3)/2., 1.),
	(3, 'l', (3,1,+1), lambda azimuth, elevation : numpy.sqrt(3./8)*numpy.cos(elevation)*(5*numpy.sin(elevation)**2 -1)*numpy.cos(azimuth), numpy.sqrt(45./32.)),
	(3, 'm', (3,1,-1), lambda azimuth, elevation : numpy.sqrt(3./8)*numpy.cos(elevation)*(5*numpy.sin(elevation)**2 -1)*numpy.sin(azimuth), numpy.sqrt(45./32.)),
	(3, 'n', (3,2,+1), lambda azimuth, elevation : numpy.sqrt(15./4)*numpy.sin(elevation)*(numpy.cos(elevation)**2)*numpy.cos(2*azimuth), numpy.sqrt(9./5.)),
	(3, 'o', (3,2,-1), lambda azimuth, elevation : numpy.sqrt(15./4)*numpy.sin(elevation)*(numpy.cos(elevation)**2)*numpy.sin(2*azimuth), numpy.sqrt(9./5.)),
	(3, 'p', (3,3,+1), lambda azimuth, elevation : numpy.sqrt(5./8)*(numpy.cos(elevation)**3)*numpy.cos(3*azimuth), numpy.sqrt(8./5.)),
	(3, 'q', (3,3,-1), lambda azimuth, elevation : numpy.sqrt(5./8)*(numpy.cos(elevation)**3)*numpy.sin(3*azimuth), numpy.sqrt(8./5.)),
] # taken from http://www.york.ac.uk/inst/mustech/3d_audio/higher_order_ambisonics.pdf
# http://ambisonics.iem.at/xchange/format/a-first-proposal-for-the-format
# http://ambisonics.iem.at/xchange/format/ambisonics-xchange-format-appendix


def nProduct(N,M) :
	import operator
	return reduce(operator.mul, xrange(N+1, N+M+1))



if __name__ == "__main__" :
	import glob
	for databaseFile in glob.glob(os.path.join(hrtf_path(),"*hrtfs")) :
		database  = HrtfDatabase(databaseFile)
		print databaseFile, database
		print databaseFile, database.layout()



