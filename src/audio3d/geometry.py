#!/usr/bin/python

import numpy as np

# this code is duplicated in scripts/AngleNormalizationTest.py,
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

# TODO: Untested
def xyz2ead(x, y, z) :
	"""Returns spherical coordinates for x,y,z in a tuple azimuth, elevation, distance"""

	d = np.sqrt(x*x+y*y+z*z)

	if (d < 0.00001):
		return 0, 0, 0

	return np.degrees(np.arctan2(y,x)), np.degrees(np.arcsin(z/d)), d


def chordDistance(e1, a1, e2, a2) :
	ra1, re1 = np.radians(a1), np.radians(e1)
	ra2, re2 = np.radians(a2), np.radians(e2)
	sa1, se1 = np.sin(ra1), np.sin(re1)
	sa2, se2 = np.sin(ra2), np.sin(re2)
	ca1, ce1 = np.cos(ra1), np.cos(re1)
	ca2, ce2 = np.cos(ra2), np.cos(re2)

	return 2 -2*(se1*se2 + ce1*ce2 * (sa1*sa2 + ca1*ca2))

	# This is a sub-optimized and less numerically accurate formula
	# left here as reference.
	dz = se1 -se2
	dy = ce1*sa1 -ce2*sa2
	dx = ce1*ca1 -ce2*ca2
	return dx*dx + dy*dy + dz*dz



# TODO: Untested
def angularDistance(a1,e1,a2,e2):
	""" Implements the haversine formula (see Wikipedia) """
	return 2*np.arcsin(np.sqrt(
			+ np.sin((e1 - e2)/2)**2
			+ np.cos(e1)*np.cos(e2)*(np.sin((a1 - a2)/2))**2
		))


# TODO: Untested
#algorithm copied from AbsoluteCoordinates2RelativeAngles.hxx (CLAM/plugins/spacialization)
def absoluteCoordinates2RelativeAngles(
			listenerX,
			listenerY,
			listenerZ,
			listenerAzimuthDegrees,
			listenerElevationDegrees,
			listenerRollDegrees,
			sourceX,
			sourceY,
			sourceZ,
		) :
	listenerAzimuth = math.radians(listenerAzimuthDegrees)
	listenerElevation = math.radians(listenerElevationDegrees)
	listenerRoll = math.radians(listenerRollDegrees)
	dx = sourceX - listenerX
	dy = sourceY - listenerY
	dz = sourceZ - listenerZ
	ca = math.cos(listenerAzimuth)
	sa = math.sin(listenerAzimuth)
	ce = math.cos(listenerElevation)
	se = math.sin(listenerElevation)
	cr = math.cos(listenerRoll)
	sr = math.sin(listenerRoll)
	seca = se*ca
	sesa = sa*se
	rotatedX = dx*ca*ce + dy*sa*ce + dz*se
	rotatedY = dx*(-sa*cr -seca*sr) + dy*(ca*cr - sesa*sr) + dz*(ce*sr)
	rotatedZ = dx*(-seca*cr +sa*sr) + dy*(-sesa*cr - ca*sr) + dz*(ce*cr)
	azimuth=math.atan2(rotatedY,rotatedX)
	module=math.sqrt(rotatedX*rotatedX+rotatedY*rotatedY+rotatedZ*rotatedZ)
	if module < 1e-12:
		elevation=rotatedZ
	else:
		elevation=math.asin(rotatedZ/module)
	return math.degrees(azimuth), math.degrees(elevation), module



