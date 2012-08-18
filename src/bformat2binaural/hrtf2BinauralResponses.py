#!/usr/bin/env python

import glob, os, re, math, numpy, sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../src/libs/python"))
import bmaudio

compensating = '--compensate' in sys.argv
verbose = '--verbose' in sys.argv or '-v' in sys.argv

higherOrder=1
if '--order' in sys.argv :
	higherOrder = int(sys.argv[sys.argv.index('--order')+1])

# A normalization factor that depends on the number of speakers
relativeNormalization = 0
if '--relative-normalization' in sys.argv :
	relativeNormalization = float(sys.argv[sys.argv.index('--relative-normalization')+1])

databaseFile = bmaudio.selectHrtfDatabase(sys.argv)

print "Using", databaseFile, "database."

sphericalHarmonics = [ #  order, name, (m,n,rho), function, weight
	(0, 'w', (0,0,+1), lambda azimuth, elevation : 1, 1./math.sqrt(2)),
	(1, 'x', (1,1,+1), lambda azimuth, elevation : math.cos(elevation)*math.cos(azimuth), 1.),
	(1, 'y', (1,1,-1), lambda azimuth, elevation : math.cos(elevation)*math.sin(azimuth), 1.),
	(1, 'z', (1,0,+1), lambda azimuth, elevation : math.sin(elevation), 1.),
	(2, 'r', (2,0,+1), lambda azimuth, elevation : (3*math.sin(elevation)**2-1)/2, 1.),
	(2, 's', (2,1,+1), lambda azimuth, elevation : math.sqrt(3./4)*math.cos(2*elevation)*math.cos(azimuth), 2./math.sqrt(3)),
	(2, 't', (2,1,-1), lambda azimuth, elevation : math.sqrt(3./4)*math.cos(2*elevation)*math.sin(azimuth), 2./math.sqrt(3)),
	(2, 'u', (2,2,+1), lambda azimuth, elevation : math.sqrt(3./4)*math.cos(elevation)**2*math.cos(2*azimuth), 2./math.sqrt(3)),
	(2, 'v', (2,2,-1), lambda azimuth, elevation : math.sqrt(3./4)*math.cos(elevation)**2*math.sin(2*azimuth), 2./math.sqrt(3)),
] # taken from http://www.york.ac.uk/inst/mustech/3d_audio/higher_order_ambisonics.pdf
# x front, y left, z up, azimuth counter-clock? starting at y?, elevation from the horizont up +-pi/2


print "Gathering files..."
hrtfDatabase = bmaudio.HrtfDatabase(databaseFile)
nImpulseResponses = len(hrtfDatabase._data)

lowerElevation = hrtfDatabase.lowerElevation()
nMissingAngles, pCompensation, vCompensation = hrtfDatabase.capCompensation(lowerElevation) if compensating else (0,1.,1.)

print nImpulseResponses, "impulse responses processed."
nCompensatedResponses = nImpulseResponses+nMissingAngles
print nCompensatedResponses, "considering compensation."

print "Processing database..."
E={}
for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
	from math import radians, degrees
	elevation = radians(elevationDegrees)
	azimuth = radians(azimuthDegrees)
	compensation = (1.,1.,1.)
	if elevation == float(lowerElevation)*math.pi/180.:
		print "c",
		compensation = (pCompensation, vCompensation,1.)
	# print "elevation:", elevation, 'azimuth:', azimuth, "x:", x, "y:", y, "z:", z, "len:", len(data)
	samplingRate, data = bmaudio.loadWave(response)
	for order, name, coefs, sphericalFunction, normFactor in sphericalHarmonics :
		if order > higherOrder : continue
		try :
			E[name] += data * sphericalFunction(azimuth,elevation) * compensation[order]
		except:
			E[name]  = data * sphericalFunction(azimuth,elevation) * compensation[order]
	sys.stdout.write(".")
	sys.stdout.flush()
print

print "Max values for"," ".join( "%s:%f"%(name,max(abs(data))) for name, data in E.iteritems())

normalizationFactor = nCompensatedResponses
if relativeNormalization :
	normalizationFactor = relativeNormalization*nImpulseResponses

for name in E.keys() :
	E[name]=E[name]/normalizationFactor

for name, data in E.iteritems() :
	if max(abs(data))<=1. :  continue
	print >> sys.stderr, "Clipping, normalization is too small. Used %f, required %f"%(
		normalizationFactor,
		normalizationFactor*max(abs(data))
		)
	sys.exit(-1)
	

for name, data in E.iteritems() :
	bmaudio.saveWave("E%s.wav"%name, data, samplingRate, verbose=True)

for name in "xyz" :
	bmaudio.saveWave("%sM.wav"%name, E['w']+E[name], samplingRate, verbose=False)
	bmaudio.saveWave("%sm.wav"%name, E['w']-E[name], samplingRate, verbose=False)


