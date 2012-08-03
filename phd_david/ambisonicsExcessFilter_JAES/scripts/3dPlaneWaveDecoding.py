#!/usr/bin/env python
from __future__ import division

import glob, os, re, math, numpy
import bmaudio
import pylab

from parameters import *
import headDistortion

from numpy import pi, sin, cos, mgrid
#from enthought.mayavi import mlab

#import pylab
#import matplotlib.axes3d
#fig = pylab.figure(figsize=(6,4))
#ax = matplotlib.axes3d.Axes3D(fig)   

import sympy

ordersToShow=[0,1,2,3,4,5,9]

def figurePath(file, extension) :
        return os.path.join('figures', os.path.splitext(os.path.basename(file))[0] + "." + extension)

def sh_normalization(l,am) :
	e = 1 if am==0 else 2
	return sympy.sqrt(sympy.Rational(e*(2*l+1)*sympy.factorial(l-am),sympy.factorial(l+am)*4)/numpy.pi) # The real one

def sh_normalization2(l,am) :
	"N_l_|a|^2 * 4pi"
	e = 1 if am==0 else 2
	return sympy.Rational(e*(2*l+1)*sympy.factorial(l-am),sympy.factorial(l+am))

def acustic_sh(l,m,a,z) :
	am = numpy.abs(m)
	if m<0 :
		return sh_normalization(l,am) * sympy.sin(am*a) * sympy.assoc_legendre(l,am,sympy.cos(z))
	else:
		return sh_normalization(l,am) * sympy.cos(am*a) * sympy.assoc_legendre(l,am,sympy.cos(z))

def projection(a0, z0, order) :
	proj = 0
	for l in xrange(0,order+1) :
		for m in xrange(-l,+l+1) :
			proj = proj + acustic_sh(l,m,a,z) * acustic_sh(l,m,a0,z0)
	return proj

def plotSpectrumsAsGreyMatrix(name, spectrums, spectralRange, yaxis, yaxisLabel="Azimuth (degrees)", levels=30) :
	pylab.rcParams['figure.figsize']=(15,7)
	pylab.rcParams['figure.subplot.left'] = 0.05
	pylab.rcParams['figure.subplot.right'] = 0.90
	spectralBinPositions = numpy.arange(0,spectralRange, spectralRange/spectrumBins)
	pylab.contourf(spectralBinPositions, yaxis, spectrums,levels,cmap=pylab.cm.Greys_r, extend='both')
	pylab.colorbar().set_label("Magnitude (dB)")
	pylab.contour(spectralBinPositions, yaxis, spectrums,levels,cmap=pylab.cm.Greys)
	pylab.xlabel("Frequency (Hz)")
	pylab.ylabel(yaxisLabel)
	pylab.savefig(name, format="pdf")
	pylab.show()

class Ambisonics3dPolarPattern() :
	def __init__(self, decoding) :
		self.x = sympy.symbols("x")
		self.pattern = sum( (decoding[l] * sh_normalization2(l,0) * sympy.legendre(l,self.x) for l in xrange(0,len(decoding)) ))
		print self.pattern
	def __call__(self, value) :
		return float(self.pattern.subs(self.x, value))


a, z = sympy.symbols("a z")
a0, z0 = 0, numpy.pi/2
epsilon = a*1e-15 + z*1e-15
order = 5
nSamples = 70
angles = numpy.arange(0,nSamples+1)*2*numpy.pi/nSamples
w = 2*math.pi*spectralRange/spectrumBins * numpy.arange(spectrumBins)

if False:
	print "Checking that the azimuthal formula is equivalent to the zenital formula"
	for l in xrange(0,order+1) :
	#	pattern1 = sum([sh_normalization2(l,m)*sympy.assoc_legendre(l,m,sympy.cos(z))*sympy.assoc_legendre(l,m,1) for m in xrange(0,l+1)])
	#	pattern1 = sh_normalization2(l,0)*sympy.assoc_legendre(l,0,sympy.cos(z))
		pattern1 = sh_normalization2(l,0)*sympy.legendre(l,sympy.cos(z))
		pylab.polar(angles, [ pattern1.subs(z,angle)/(2*l+1) for angle in angles ], label="Zenital %s"%l)

		pattern2 = sum([sh_normalization2(l,m)*sympy.cos(m*z)*sympy.assoc_legendre(l,m,0)**2 for m in xrange(0,l+1)])
		pylab.polar(angles, [ pattern2.subs(z,angle)/(2*l+1) for angle in angles ], label="Azimuthal %s"%l)
	pylab.title("Zenital vs Azimuthal variation",horizontalalignment='center', verticalalignment='baseline', position=(.5,-.1))
	pylab.rgrids(numpy.arange(.4,1,.2),angle=220)
	pylab.legend(loc=2)
	pylab.savefig(figurePath(__file__,"pdf"))
	pylab.show()
# We take the azimuthal formula which is faster
print "Computing component patterns"
patternComponents = [
	sh_normalization2(l,0)*sympy.legendre(l,sympy.cos(z))
	for l in xrange(0,order+1)
	]
for l, pattern in enumerate(patternComponents) :
	print "%i:"%l, pattern

if False:
	print "Combining components to get the decoding pattern"
	for l in xrange(0,order+1) :
		decoding = headDistortion.decoding3dInPhase(l)/(l+1)
		decoding = headDistortion.decoding3dBasic(l)/(l+1)**2
		decoding = headDistortion.decoding3dMaxRe(l)/math.sqrt(1./(l+1))
		print l, decoding
#		pattern = sum([g*component for g,component in zip(decoding,patternComponents[:l+1])])
#		patternValues = numpy.array([ float(pattern.subs(z,angle)) for angle in angles ])
		pattern = Ambisonics3dPolarPattern(decoding)
		patternValues = numpy.array( [ pattern(math.cos(angle)) for angle in angles ] )
		print "max", max(patternValues)
		pylab.polar(angles, patternValues, label="%s"%l)
	pylab.title("InPhase",horizontalalignment='center', verticalalignment='baseline', position=(.5,-.1))
	pylab.rgrids(numpy.arange(.4,1,.2),angle=220)
	pylab.legend(loc=2)
	pylab.savefig(figurePath(__file__,"pdf"))
	pylab.show()

if False:
	print "Filtering at 90 degrees"
	plot = bmaudio.SpectrumDisplay()
	plot.inDb()
	#plot.showPhase()
	plot.ylim(-40,5)
	delays = [
		(bmaudio.sphericalHeadDelay(azimuth, R, r)-(r-R))/c
		for azimuth in angles
		]
	for no in ordersToShow[:-1] :
		print no
		decoding = headDistortion.decoding3dMaxRe(no) / (no+1)
		decoding = headDistortion.decoding3dInPhase(no) 
		amplitudes = [ sum(
				(decoding[l] * float(sh_normalization2(l,0)) * float(sympy.legendre(l,cosTheta)) * math.sin(azimuth) / nSamples
				for l in xrange(0,no+1) ) )
			for azimuth, cosTheta in zip(angles, numpy.cos(numpy.pi/2-angles)) ]
		plot.addSpectrumData(headDistortion.delayAmplitudeToSpectrum(delays,amplitudes,w),spectralRange,"Order %i"%no)
	plot.show()

def deviated3dAmbisonics(pattern, planeWaveAzimuth, ringAzimuth, slices=72) :
	"""
	Discretely integrates a ring perpendicular to the Y axis of the decoding of a plane wave
	planeWaveAzimuth: angle that forms the incoming wave direction and the left ear (azimuth 90)
	ringAzimuth: azimuth at which we are integrating the ring
	"""
	return sum(
		( float(pattern(cosTheta))
			for cosTheta in 
				+math.sin(ringAzimuth)*math.sin(planeWaveAzimuth)
				-math.cos(ringAzimuth)*math.cos(planeWaveAzimuth)*numpy.cos(numpy.arange(0,2*pi,2*pi/slices))
			)) / slices  # TODO: I think that that should be cos but i get zeros.


if True : 
	no = 3
	degreesResolution=10
	# the angle between the incomming plane wave and the interaural axis
	decoding = headDistortion.decoding3dMaxRe(no)
	decoding = headDistortion.decoding3dInPhase(no)
	halfAngles = numpy.linspace(-90, 90, num=nSamples)
	planeWaveAzimuths = xrange(-90,90+1,degreesResolution)
	spectrums = numpy.zeros((len(planeWaveAzimuths), spectrumBins))
	delays = [
		(bmaudio.sphericalHeadDelay(azimuth, R, r)-(r-R))/c
		for azimuth in halfAngles
		]
	pattern = Ambisonics3dPolarPattern(decoding)
	for i,pwAzimuth in enumerate(planeWaveAzimuths) :
		print pwAzimuth
		magnitudes = [ deviated3dAmbisonics(pattern, math.radians(pwAzimuth), azimuth, slices=360)/len(halfAngles) for azimuth in halfAngles ]
		spectrum = headDistortion.delayAmplitudeToSpectrum(delays, magnitudes, w)
		spectrums[i,:] = 20*numpy.log(numpy.abs(spectrum))
	print spectrum.size, planeWaveAzimuths

	plotSpectrumsAsGreyMatrix("boo.pdf", spectrums, spectralRange, planeWaveAzimuths, yaxisLabel="Azimuth (degrees)")#, levels=numpy.arange(-60,6+1,6))


if False :
	print "Component filter"
	plot = bmaudio.SpectrumDisplay()
	plot.inDb()
	plot.showPhase()
	#plot.ylim(-40,5)
	for no in ordersToShow[:-1] :
		print no
		spectrum = numpy.zeros(spectrumBins, numpy.complex)
		for azimuth in angles :
			t = (bmaudio.sphericalHeadDelay(azimuth, R, r)-(r-R))/c
			cosTheta = math.cos(math.pi/2-azimuth)
			amplitude = float(sh_normalization2(no,0)) * float(sympy.legendre(no,cosTheta)) * math.cos(azimuth) / nSamples
			sinusoid = numpy.exp( 1j * t * w )
			spectrum += amplitude * sinusoid
		plot.addSpectrumData(spectrum,spectralRange,"Order %i"%no)
	plot.show()


if False :
	print "3D representation of the plane wave decoding spherical pattern"
	p=sympy.Plot()
	p.append(epsilon+numpy.abs(projection(2,0*numpy.pi/4,order)), "mode=spherical;color=zfade4")
	p.append(epsilon+numpy.abs(projection(2,1*numpy.pi/4,order)), "mode=spherical;color=rainbow")
	p.append(epsilon+numpy.abs(projection(2,2*numpy.pi/4,order)), "mode=spherical;color=zfade")
	p.append(epsilon+numpy.abs(projection(2,3*numpy.pi/4,order)), "mode=spherical;color=zfade3")
#	p.append(numpy.abs(sum([(sh_normalization(o,l)**2 for l,o in (xrange(0,o+1)) for o in xrange(0,l+1) ])))
#	p.append(epsilon+(projection(2+2*numpy.pi/4,3*numpy.pi/4,order)), "mode=spherical;color=zfade3")
#	p.append(epsilon+(projection(2+3*numpy.pi/4,3*numpy.pi/4,order)), "mode=spherical;color=zfade3")
#	p.append(epsilon+40*sympy.abs(acustic_sh(1,-1,a,z)), "mode=spherical;color=zfade4")
	p.show()
	input

print "Plm(0) Plm(1)"
for l in xrange(0,order+1) :
	for m in xrange(0,l+1) :
		print l,m,":",sympy.assoc_legendre(l,m,0), sympy.assoc_legendre(l,m,1)







