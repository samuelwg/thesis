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




if __name__ == "__main__" :
	sys.exit(unittest.main())

