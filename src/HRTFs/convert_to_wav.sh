#!/bin/bash
unzip KEMAR_HRTF.zip
for dir in elev*
do 
	for a in $dir/*.dat
	do 
		wave=$dir/`basename $a .dat`.wav
		sox -t raw -r44100 -w -s -x $a $wave
	done
done
