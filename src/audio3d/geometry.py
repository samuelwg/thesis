#!/usr/bin/python



def normalizeAngles(azimut, elevation) :
	elevation += 90.
	elevation %= 360.
	if elevation > 180. :
		azimut += 180.
		elevation = 360 - elevation
	elevation -= 90.
	azimut %= 360.
	return azimut,elevation



