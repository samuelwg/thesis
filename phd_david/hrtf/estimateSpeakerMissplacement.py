#!/usr/bin/env python

import glob, os, re, math, numpy, sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../../../src/libs/python"))
from bmaudio import *

databaseFile = selectHrtfDatabase(sys.argv)
print "Using", databaseFile, "database."
print "Gathering files..."
hrtfDatabase = HrtfDatabase(databaseFile)

deviations=[]
sampleDelays=[]


for elevationDegrees, azimuthDegrees, response in hrtfDatabase._data :
	from math import radians, degrees
	elevation = radians(elevationDegrees)
	azimuth = radians(azimuthDegrees)
	z = math.sin(elevation)
	y = math.cos(elevation) * math.sin(azimuth)
	x = math.cos(elevation) * math.cos(azimuth)
	samplingRate, data = loadWave(response)
	sampleDelay = firstMaxima(data,.4)
	sampleDelays.append((x,y,z,sampleDelay, azimuth, elevation))

removedSamples = 140
# TODO: This is to compute the speaker displacement error distribution
# from http://www.grasinfo.dk/documents/KEMAR_Manikin_Measurements.PDF
# Interaural distance: 152mm, 
earDistance = .152/2
aboveCutoff=40
belowCutoff=-40

def centralizeDelay(maximaLocation, x, y, z, h) :
	return math.sqrt(maximaLocation**2 - h**2*(x**2+z**2)) - h*y

directDelays = [centralizeDelay((sampleDelay+removedSamples)*c/samplingRate, x, y, z, earDistance)
	for x,y,z,sampleDelay, azimuth, elevation in sampleDelays 
		if elevation<math.radians(aboveCutoff) and elevation>math.radians(belowCutoff) and
		azimuth>math.radians(190) and azimuth<math.radians(260) ]

for x,y,z,sampleDelay, azimuth, elevation in sampleDelays  :
	if elevation<math.radians(aboveCutoff) and elevation>math.radians(belowCutoff) and \
	azimuth>math.radians(190) and azimuth<math.radians(260) :
		if centralizeDelay((sampleDelay+removedSamples)*c/samplingRate, x, y, z, earDistance) < 1.46:
			print math.degrees(azimuth), math.degrees(elevation)
		if centralizeDelay((sampleDelay+removedSamples)*c/samplingRate, x, y, z, earDistance) > 1.50:
			print math.degrees(azimuth), math.degrees(elevation)

import stats
print "Direct delays stdev:", stats.stdev(directDelays)
print "Direct delays mean:", stats.mean(directDelays)

import pylab
pylab.grid()
pylab.hist(directDelays,bins=7)
pylab.show()

nominalDistance = 1.4
meanDelay = sum((sampleDelay for x,y,z,sampleDelay,_,_ in sampleDelays))*c/samplingRate/len(sampleDelays)
for x,y,z,sampleDelay,azimuth,elevation in sampleDelays :
	distanceDelay = sampleDelay * c / samplingRate - meanDelay + nominalDistance
	centralizedDelay = centralizeDelay(distanceDelay, x, y, z, earDistance)
	deviations.append(sampleDelay*c/samplingRate)

	if (abs(sampleDelay-41)>13) : 
		print response, sampleDelay
	if (False and abs(centralizedDelay-1.4)>0.03) : 
		print response, sampleDelay, distanceDelay, centralizedDelay, centralizedDelay-nominalDistance


	# print "elevation:", elevation, 'azimuth:', azimuth, "x:", x, "y:", y, "z:", z, "len:", len(data)
	sys.stdout.write(".")
	sys.stdout.flush()

deviations = numpy.array(deviations)
print
print "Speaker distance deviations. Max:", deviations.max(), "min:", deviations.min(), "std:", deviations.std(), "mean:", deviations.mean()



