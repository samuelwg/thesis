#!/usr/bin/python

# Install wav2png:
#  $ sudo apt-get install libgd2-xpm-dev
#  svn checkout http://wav2png.googlecode.com/svn/trunk/ wav2png
#  cd build/linux
#  $ sudo cp wav2png /usr/local/bin/ 

import os, sys, string
sys.path.append('../../clam/CLAM/scripts/')
from audiob2b import runBack2BackProgram

data_path='b2b/bformat2binaural'
if not os.access(data_path, os.X_OK) :
	print "\033[31mData path not found:\033[0m ", data_path
	print """\
You need a to checkout the data repository and link its b2b/ dir. 
For example:
$ cd
$ svn co http://parumi.org/acustica/data data_acustica
$ cd acustica/visualitzador_escena_c++
$ ln -s ~/data_acustica/b2b/
"""
	sys.exit(-1) 

back2backTests = [
	(prefix, commandLine, ['E%s.wav'%component for component in 'wxyz'] )
	for prefix, commandLine in [
		('MIT_noncompensated',
			'./hrtf2BinauralResponses.py '),
		('MIT_noncompensated_R',
			'./hrtf2BinauralResponses.py -r'),
		('MIT_compensated',
			'./hrtf2BinauralResponses.py --compensate'),
		('MITdiffuse_noncompensated',
			'./hrtf2BinauralResponses.py --mitdiffuse'),
		('MITdiffuse_compensated',
			'./hrtf2BinauralResponses.py --mitdiffuse --compensate'),
		('IRCAM1002_noncompensated',
			'./hrtf2BinauralResponses.py --ircam 1002'),
		('IRCAM1002_compensated',
			'./hrtf2BinauralResponses.py --ircam 1002 --compensate'),
		('IRCAM1002R_noncompensated',
			'./hrtf2BinauralResponses.py --ircam 1002 -r'),
		('IRCAM1002R_compensated',
			'./hrtf2BinauralResponses.py --ircam 1002 -r --compensate'),
		('KREUZER_noncompensated',
			'./hrtf2BinauralResponses.py --kreuzer'),
		('KREUZER_compensated',
			'./hrtf2BinauralResponses.py --kreuzer --compensate'),
	]]

runBack2BackProgram(data_path, sys.argv, back2backTests)




