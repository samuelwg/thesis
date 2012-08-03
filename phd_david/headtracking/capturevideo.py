#!/usr/bin/python


import pygame
import Image
from pygame.locals import QUIT, KEYDOWN
import sys

import opencv
import cv
#this is important for capturing/displaying images
from opencv import highgui 
import ImageChops
import ImageOps
import ImageDraw

import freenect
import numpy


cascadeFile = "/usr/share/doc/opencv-doc/examples/haarcascades/haarcascades/haarcascade_frontalface_default.xml.gz"
cascadeFile = "/usr/share/doc/opencv-doc/examples/haarcascades/haarcascades/haarcascade_frontalface_alt.xml.gz"
cascadeFile = "/usr/share/doc/opencv-doc/examples/haarcascades/haarcascades/haarcascade_frontalface_alt2.xml.gz"
cascadeFile = "/usr/share/doc/opencv-doc/examples/haarcascades/haarcascades/haarcascade_eye_tree_eyeglasses.xml.gz"
#cascade = opencv.cv.cvLoad(cascadeFile)
cascade = opencv.cv.cvLoadHaarClassifierCascade(cascadeFile, opencv.cv.cvSize(1,1))
camera = highgui.cvCreateCameraCapture(0)
storage = opencv.cvCreateMemStorage(0)

fps = 10.0
size = (640,480)
pygame.init()
window = pygame.display.set_mode(size)
pygame.display.set_caption("WebCam Demo")
screen = pygame.display.get_surface()
black = Image.new('RGB', size)
red = Image.new('RGB', size, (255,0,0))

def movementDetector(prev, current) :
	# meassure the difference
	step = ImageChops.difference(prev, current)
	# join channels
	step = step.convert("L")
	# thresholding
	step = Image.eval(step, lambda p: 100 if p>10 else 0)
	# compositing with current
	step = Image.composite(red, current, step)
	return step

def detectFaces(image) :
	return []
	opencv.cv.cvClearMemStorage(storage)
	image_scale=4.0
	gray = opencv.cv.cvCreateImage( opencv.cv.cvSize(image.width,image.height), 8, 1 )
	small_img = opencv.cv.cvCreateImage(
		(
			opencv.cv.cvRound(image.width/image_scale),
			opencv.cv.cvRound(image.height/image_scale)
		),
		8, 1 )

	# convert color input image to grayscale
	opencv.cv.cvCvtColor( image, gray, opencv.cv.CV_BGR2GRAY )

	# scale input image for faster processing
	opencv.cv.cvResize( gray, small_img, opencv.cv.CV_INTER_LINEAR )
	opencv.cv.cvEqualizeHist( small_img, small_img )
	faces = opencv.cv.cvHaarDetectObjects(
		small_img,
		cascade,
		storage,
		1.1,
		2,
		opencv.CV_HAAR_DO_CANNY_PRUNING,
		opencv.cv.cvSize(5, 5),
#		opencv.cv.cvSize(400, 400),
		)
	return [(
		opencv.cv.cvRound(box.x*image_scale),
		opencv.cv.cvRound(box.y*image_scale),
		opencv.cv.cvRound((box.x+box.width)*image_scale),
		opencv.cv.cvRound((box.y+box.height)*image_scale),
		) for box in faces]

def addBoxes(image, boxes) :
	draw = ImageDraw.Draw(image)
	for box in boxes:
		print "Detected:", box
		draw.ellipse(box, outline=(0,0,100))
 

prev=None

print "Action!"
while True:
	for event in pygame.event.get() :
		if event.type == QUIT : sys.exit(0)
	cvimage = highgui.cvQueryFrame(camera)
	image = opencv.adaptors.Ipl2PIL(cvimage)
	depthImage = freenect.sync_get_depth()
	image = cv.CreateImage(
		depthImage.shape[1],depthImage.shape[0],
		cv.IPL_DEPTH_16U, depthImage.shape[2],
		depthImage)
	faces = detectFaces(cvimage)
	if prev is None : prev = image
	result = movementDetector(prev, image)
	addBoxes(result, faces)
	pg_img = pygame.image.frombuffer(result.tostring(), result.size, result.mode)
	screen.blit(pg_img, (0,0))
	pygame.display.flip()
#	pygame.time.delay(int(1000 * 1.0/fps))
	prev = image



