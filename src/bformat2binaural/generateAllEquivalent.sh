#!/bin/bash

# This script generate equivalents for all supported databases and places them
# in a different folder for each database.

function die()
{
	echo $1 >&2
	exit $2
}

function moveEquivalentsToDir()
{
	echo "Moving to $1"
	mkdir -p $1
	mv E{r,s,t,u,v,w,x,y,z}.wav $1
}


for database in ../HRTFs/*.hrtfs
do
	echo Processing $database
	./hrtf2BinauralResponses.py --order 2 --compensate --database $database || die "Error generating equivalents for database $database" -1
	moveEquivalentsToDir equivalentIRs/`basename $database .hrtfs`
done


