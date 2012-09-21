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

"""
TODO:
"""

import os, sys
from geometry import normalizeAngles

def hrtf_path() :
	if 'HRTF_PATH' not in os.environ.keys():
		print >> sys.stderr, 'Could not find HRTF_PATH environment variable. Please define it and point it to where imm_bm/HRTFs is'
		sys.exit(-1)
	return os.environ['HRTF_PATH']

def selectHrtfDatabase( args ) :

	if '--help' in args :
		print >> sys.stderr, """
	--database PATH     Explicitly choses a database path
	-r                  Choose the right ear (default left)
	--speakers N        Choose a 2d sublayout of N speakers
	--mitdiffuse        Choose the difuse filtered MIT database (default Full MIT)
	--mitcompact        Choose the compact MIT database (default Full MIT)
	--ircam N           Choose the IRCAM database for subject N (default Full MIT)
    --kreuzer           Choose the Kreuzer simulated database (default (Full MIT)
"""

	if '--database' in args :
		return args[args.index('--database')+1]

	# choose ear, left default
	earLetter = 'R' if '-r' in args else 'L'

	speakerInfix = ''
	if '--speakers' in args :
		# TODO: Interleaved should work again
		density = int(args[args.index('--speakers')+1]) # TODO: Remove it
		speakerInfix="-2d%02i"%density

	# IRCAM database place left at 90 degrees, no sublayouts
	if '--ircam' in args :
		nIrcam = int(args[args.index('--ircam')+1])
		return os.path.join(hrtf_path(),"ircam%i%s.hrtfs"%(nIrcam, earLetter))

	# KREUZER database place left at 90 degrees, no left/right
	if '--kreuzer' in args :
		return os.path.join(hrtf_path(),"kreuzerDatabase2%s.hrtfs"%(speakerInfix))

	# MIT KEMAR database (full, compact and difuse)
	if '--mitdiffuse' in args :
		return os.path.join(hrtf_path(),"mitKemarDiffuse%s%s.hrtfs" % (earLetter, speakerInfix))
	if '--mitcompact' in args :
		return os.path.join(hrtf_path(),"mitKemarCompact%s%s.hrtfs" % (earLetter, speakerInfix))
	return os.path.join(hrtf_path(),"mitKemarFull%s%s.hrtfs" % (earLetter, speakerInfix))


class HrtfDatabase(object) :

	def __init__(self, databaseFile) :
		self._databaseFile = databaseFile
		self._data = []
		base = os.path.dirname(databaseFile)
		for line in open(databaseFile) :
			try : elevation, azimuth, filename = line.split()
			except: continue
			azimuth, elevation = normalizeAngles(float(azimuth), float(elevation))
			self._data.append( ( float(elevation), float(azimuth), os.path.join(base,filename)) )
		self._orientationToFilename = dict(((e,a),f) for e,a,f in self._data)


class HrtfDatabase2 :

	def __init__(self, databaseFile) :
		self._databaseFile = databaseFile
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

