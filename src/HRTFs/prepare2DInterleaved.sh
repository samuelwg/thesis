#!/bin/bash

# This script prepares interleaved 2D versions of
# mit (diffuse, full) and kreuzer2 databases
# for different number of speakers

function configurations
{
	echo 2 3 4 6 8 9 12 18 24 36 72
}

echo -e
for configuration in `configurations`
do
	paddedConfiguration=`printf '%02i' $configuration`
	for hrtfFile in \
		"kreuzerDatabase2-2d$paddedConfiguration" \
		"mitKemarDiffuseL-2d$paddedConfiguration"  \
		"mitKemarFullL-2d$paddedConfiguration"  \
		"mitKemarFullR-2d$paddedConfiguration"  ;
	do
		rm -f $hrtfFile.hrtfs
		touch $hrtfFile.hrtfs
	done
	gap=$((360/configuration))
	for ((angle=0; angle<360; angle+=gap))
	do
		paddedAngle=`printf '%03i' $((angle))`
		mitAngle=$((360-angle))
		echo 0 $angle kreuzer2/output_e+00_a$paddedAngle.0.wav >> kreuzerDatabase2-2d$paddedConfiguration.hrtfs
		echo 0 $mitAngle MIT_KEMAR/full/elev0/L0e$paddedAngle'a'.wav >> mitKemarFullL-2d$paddedConfiguration.hrtfs
		echo 0 $angle MIT_KEMAR/full/elev0/R0e$paddedAngle'a'.wav >> mitKemarFullR-2d$paddedConfiguration.hrtfs
		echo 0 $angle MIT_KEMAR/diffuse/elev0/L0e$paddedAngle'a'.wav >> mitKemarDiffuseL-2d$paddedConfiguration.hrtfs
	done
done

