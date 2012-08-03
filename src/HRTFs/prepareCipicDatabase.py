#!/usr/bin/python
import glob
import scipy.io
import re
import numpy
import os
import sys
sys.path.append('../src/libs/python/')
import bmaudio

def run(command):
	print "\033[34m"+command+"\033[0m"
	error = os.system(command)
	if error :
		print "\033[31mFailed!\033[0m"
		sys.exit(error)
	print "\033[33mOK!\033[0m"
def norun(command) : 
	pass

run("wget -c 'http://interface.cipic.ucdavis.edu/data/CIPIC_hrtf_database.zip'")
run("unzip CIPIC_hrtf_database.zip")

files = glob.glob("CIPIC_hrtf_database/standard_hrir_database/subject_*/hrir_final.mat")

elevations=-45+5.625*numpy.array(range(0,50))
azimuths=[-80, -65, -55]+range(-45,46,5)+[55, 65, 80]
print len(azimuths)
print elevations
print azimuths


for matfile in files :
	matdata = scipy.io.loadmat(matfile)
	subject = re.findall('[0-9]+',matdata['name'][0])[0] # has the form subject_000
	os.system('mkdir -p "CIPIC/subject_%s"' % subject)
	for channel in ['l','r'] :
		database=open("cipic%s%s.hrtfs"%(subject,channel.upper()),'w')
		hrirs = matdata['hrir_'+channel]
		normFactor = max(abs(hrirs.reshape(hrirs.size)))
		hrirs /= normFactor
		for azimuthIndex, azimuthData in enumerate(hrirs) :
			for elevationIndex, hrir in enumerate(azimuthData) :
				azimuth = azimuths[azimuthIndex]
				elevation = elevations[elevationIndex]
				print "Subject", subject, "original", elevation, azimuth, "->",
				if elevation > 90 :
					azimuth = 180 - azimuth
					elevation = 180 - elevation
				if channel == 'l' : azimuth = (360 - azimuth)%360
				print elevation, azimuth
				wavefile = 'CIPIC/subject_%s/%s_e%+03.1f_a%+03i.wav'%(subject,channel.upper(),elevation,azimuth)
				print >> database, elevation, azimuth, wavefile
				bmaudio.saveWave(wavefile, hrir, 44100)



