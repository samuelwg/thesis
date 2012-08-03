import numpy

N=72
r=1.4
R=0.075 # interaural distance
R=0.088 # equivalent sphere
nBins = 1024*8
samplingRate = 44100
spectralRange = samplingRate/2.
c = 344.
spectrumBins = nBins/2+1
ordersToShow = numpy.arange(0,100,20)
ordersToShow = numpy.array([0, 1, 2, 3, 5, 9, 35])
#ordersToShow = numpy.arange(13)*3+36
#ordersToShow = numpy.array([7, 41, 51, 61, 68, 69, 70, 71, 72])
#ordersToShow = 2*numpy.array([7, 41, 51, 61, 68, 69, 70, 71, 72])
#ordersToShow = numpy.array([N-3, N-2, N-1, N, N+1]) # to see what happens with N multiples -1
NOrder = max(ordersToShow)+1



