#!/usr/bin/python
"""
This scripts can be used to find speaker configuration with the following criterium:
- Speakers are layed out in horizontal rings at monotone elevations,
  being one of the rings in the azimuth 0
- For each ring, as many speakers are placed so that it has the same
  speaker density (speaker/surface) than the azimuth 0 ring
Given an elevation step in degrees it tries all the even number of speakers 
at the zero azimuth and tries to match that density on all the levels.
An error on the density missmatch is given for each configuration.
"""


import math
import numpy
import sys

elevationStep = int(sys.argv[1])
if 180%elevationStep :
	print >> sys.stderr, "The elevation step should be a divisor of 180"
	sys.exit()
if elevationStep%2 :
	print >> sys.stderr, "This elevation step does not render a symmetrical distribution"
	sys.exit()

def surfaceForElevation(elevationStep, elevation) :
	"computes the sphere surface comprised +-elevationStep around the given elevation"
	upperIntegral = math.sin(math.radians(elevation+elevationStep/2.))
	lowerIntegral = math.sin(math.radians(elevation-elevationStep/2.))
	return 2*math.pi*(upperIntegral-lowerIntegral)

elevations = range(-90+elevationStep/2.,90,elevationStep)
elevationSurfaces = [ surfaceForElevation(elevationStep, e) for e in elevations ]

print zip(elevations, elevationSurfaces)

minSurface = min(elevationSurfaces)
oldN = []

for cte in numpy.arange(2,70,.001) :
	N = [ int(round ( surface/minSurface * cte )) for surface in elevationSurfaces ]
	if N == oldN : continue
	oldN = N
	nSpeakers = sum(N)
	speakerSurfaces = [ surface/n for surface, n in zip(elevationSurfaces, N) ]
	error = sum([ n * ( (surface/n)-1./nSpeakers)**2 for surface, n in zip(elevationSurfaces, N)]) 
	degreesDeviation = math.sqrt(sum([ n*((round(36000./n)-36000./n)**2) for n in N ]))

#	if degreesDeviation >= 1 : continue
	print "Points:", int(sum(N))
	print "degreesDeviation:" , degreesDeviation
	print "Deviation:", error
	
	for elevation, speakers in reversed(zip(elevations, N )) :
		if elevation < 0: continue
		print elevation, speakers, 360./speakers
	print


