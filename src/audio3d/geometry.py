#!/usr/bin/python

import numpy as np

def normalizeAngles(azimut, elevation) :
	elevation += 90.
	elevation %= 360.
	if elevation > 180. :
		azimut += 180.
		elevation = 360 - elevation
	elevation -= 90.
	azimut %= 360.
	return azimut,elevation

def ead2xyz(e,a,d) :
	"""
		Converts polars to cartesian with those conventions:
		x front, y left, z up
		0 elevation and azimuth front
		positive elevation up, positive azimuth left
		Angles are in degrees.
	"""
	ra, re = np.radians(a), np.radians(e)
	sa, se = np.sin(ra), np.sin(re)
	ca, ce = np.cos(ra), np.cos(re)
	x,y,z = d*ce*ca, d*ce*sa, d*se
	return x,y,z



