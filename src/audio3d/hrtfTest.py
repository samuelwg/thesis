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

from hrtf import selectHrtfDatabase, hrtf_path, HrtfDatabase
import unittest
from numpy.testing import assert_equal as np_assert_equal

import os, sys

class SelecHrtfDatabaseTest(unittest.TestCase) :

	def setUp(self) :
		try :
			self._oldHrtfPath = os.environ['HRTF_PATH']
		except KeyError :
			self._oldHrtfPath = None

	def tearDown(self) :
		if self._oldHrtfPath is None :
			del os.environ['HRTF_PATH']
		else :
			os.environ['HRTF_PATH'] = self._oldHrtfPath

	def test_hrtf_path(self) :
		os.environ['HRTF_PATH'] = 'mypath'
		self.assertEqual('mypath', hrtf_path())

	def test_hrtf_path_whenNoDefined(self) :
		if 'HRTF_PATH' in os.environ : del os.environ['HRTF_PATH']
		try :
			hrtf_path()
			self.fail("Exit expected")
		except SystemExit :
			# TODO: Check and mute the stderr message
			pass

	def test_selectHrtfDatabase(self) :
		os.environ['HRTF_PATH'] = 'basepath'
		self.assertEqual(
			'basepath/mitKemarFullL.hrtfs',
			selectHrtfDatabase("".split()))

	def test_selectHrtfDatabase_withOtherPath(self) :
		os.environ['HRTF_PATH'] = 'otherpath'
		self.assertEqual(
			'otherpath/mitKemarFullL.hrtfs',
			selectHrtfDatabase("".split()))

	def test_selectHrtfDatabase_trailingSlash(self) :
		os.environ['HRTF_PATH'] = 'basepath/'
		self.assertEqual(
			'basepath/mitKemarFullL.hrtfs',
			selectHrtfDatabase("".split()))

	def test_selectHrtfDatabase_explicit(self) :
		os.environ['HRTF_PATH'] = 'basepath'
		self.assertEqual(
			'my/database.hrtfs',
			selectHrtfDatabase("--database my/database.hrtfs".split()))

	def test_selectHrtfDatabase_right(self) :
		os.environ['HRTF_PATH'] = 'basepath'
		self.assertEqual(
			'basepath/mitKemarFullR.hrtfs',
			selectHrtfDatabase("-r".split()))

	def test_selectHrtfDatabase_2dSublayout(self) :
		os.environ['HRTF_PATH'] = 'basepath'
		self.assertEqual(
			'basepath/mitKemarFullL-2d20.hrtfs',
			selectHrtfDatabase("--speakers 20".split()))

	def test_selectHrtfDatabase_ircam(self) :
		os.environ['HRTF_PATH'] = 'basepath/'
		self.assertEqual(
			'basepath/ircam15L.hrtfs',
			selectHrtfDatabase("--ircam 15".split()))

	def test_selectHrtfDatabase_kreuzer(self) :
		os.environ['HRTF_PATH'] = 'basepath/'
		self.assertEqual(
			'basepath/kreuzerDatabase2.hrtfs',
			selectHrtfDatabase("--kreuzer".split()))

	def test_selectHrtfDatabase_mitDiffuse(self) :
		os.environ['HRTF_PATH'] = 'basepath/'
		self.assertEqual(
			'basepath/mitKemarDiffuseL.hrtfs',
			selectHrtfDatabase("--mitdiffuse".split()))

	def test_selectHrtfDatabase_mitCompact(self) :
		os.environ['HRTF_PATH'] = 'basepath/'
		self.assertEqual(
			'basepath/mitKemarCompactL.hrtfs',
			selectHrtfDatabase("--mitcompact".split()))

import shutil
import wavefile

class HrtfDatabaseTest(unittest.TestCase) :
	def setUp(self) :
		try: 
			shutil.rmtree("testhrtf")
		except OSError : pass
		os.mkdir("testhrtf")
		self.writefile("testhrtf/db.hrtfs",
			"0  0    front.wav\n"
			"0  90   left.wav\n"
			"0  270  right.wav\n"
			"0  180  back.wav\n"
			)
		import numpy as np
		def saveWave(name, data) :
			redata = data.astype(np.float32)[:,np.newaxis]
			with wavefile.WaveWriter('testhrtf/%s.wav'%name) as writer :
				writer.write(redata)
			return redata
			
		self._audioFront = saveWave('front', np.arange(.1,.3,.01))
		self._audioBack  = saveWave('back',  np.arange(.2,.4,.01))
		self._audioLeft  = saveWave('left',  np.arange(.3,.5,.01))
		self._audioRight = saveWave('right', np.arange(.4,.6,.01))

	def tearDown(self) :
		shutil.rmtree("testhrtf")

	def writefile(self, name, content) :
		f = open(name, 'w')
		f.write(content)
		f.close()

	def test_repr(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		self.assertEqual(
			"HrtfDatabase:["
				"(0.0, 0.0, 'testhrtf/front.wav'), "
				"(0.0, 90.0, 'testhrtf/left.wav'), "
				"(0.0, 270.0, 'testhrtf/right.wav'), "
				"(0.0, 180.0, 'testhrtf/back.wav')]"
			, repr(db))

	def test_existing(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		self.assertEqual('testhrtf/left.wav', db.wavefile(0,90))

	def test_missing(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		try :
			db.wavefile(0,91)
			self.fail("An exception was expected")
		except KeyError, e:
			self.assertEqual(e.message, (0,91))

	def test_nearest_exact(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		self.assertEqual((0,90), db.nearest(0, 90))

	def test_nearest_near(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		self.assertEqual((0,90), db.nearest(0, 91))

	def test_nearestWavefile(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		self.assertEqual(
			'testhrtf/left.wav',
			db.nearestWavefile(0, 91))

	def test_loadAudio_left(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		np_assert_equal(
			self._audioLeft,
			db.hrtf(0,90))

	def test_loadAudio_nearleft(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		np_assert_equal(
			self._audioLeft,
			db.hrtf(0,91))

	def test_cachedFiles_single(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		db.hrtf(0,90)
		self.assertEqual([
			"testhrtf/left.wav",
			], db.cachedFiles())

	def test_cachedFiles_many(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		db.hrtf(0,+90)
		db.hrtf(0,-90)
		self.assertEqual([
			"testhrtf/left.wav",
			"testhrtf/right.wav",
			], db.cachedFiles())

	def test_cachedFiles_same(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		db.hrtf(0,+90)
		db.hrtf(0,+91)
		self.assertEqual([
			"testhrtf/left.wav",
			], db.cachedFiles())

	def test_cachedFiles_all(self) :
		db = HrtfDatabase('testhrtf/db.hrtfs')
		db.preload()
		self.assertEqual([
			"testhrtf/left.wav",
			"testhrtf/front.wav",
			"testhrtf/back.wav",
			"testhrtf/right.wav",
			], db.cachedFiles())


if __name__ == "__main__" :
	sys.exit(unittest.main())

