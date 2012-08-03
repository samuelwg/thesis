#!/usr/bin/python

import yaml
import numpy

def opencv_matrix_constructor(loader, node) :
	value = loader.construct_mapping(node, deep=True)
	typeMapping = dict(
		d=numpy.float,
		i=numpy.int,
		)
	data = value['data']
	shape = value['rows'], value['cols']
	type = dtype=typeMapping[value['dt']]
	result = numpy.array(data, dtype=type)
	return result.reshape(shape)

yaml.add_constructor('tag:yaml.org,2002:opencv-matrix', opencv_matrix_constructor)

calibrationInfo = yaml.load(file("kinect_calibration.yml"))


for name, matrix in calibrationInfo.iteritems() :
	print name, "=\n", matrix


