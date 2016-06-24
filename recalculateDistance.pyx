import numpy
cimport numpy
cimport cython

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def recalculateDistanceFunction(numpy.ndarray[dtype=DTYPE_t, ndim=1] imageCenter,
                                numpy.ndarray[dtype=DTYPE_t, ndim=2] predictedPattern, 
                                numpy.ndarray[dtype=DTYPE_t, ndim=2] detectedPeaks,
                                float pixelSize):
                            
    cdef int detectedPeakRowIdx
    cdef int detectedPeaksNrows = detectedPeaks.shape[0]
    cdef int predictedSpotRowIdx
    cdef int predictedPatternNrows = predictedPattern.shape[0]
    cdef float calcX
    cdef float calcY
       
    for detectedPeakRowIdx in range(0, detectedPeaksNrows):       
        for predictedSpotRowIdx in range(0, predictedPatternNrows):
            if predictedPattern[predictedSpotRowIdx, 0] == detectedPeaks[detectedPeakRowIdx, 0] and predictedPattern[predictedSpotRowIdx, 1] == detectedPeaks[detectedPeakRowIdx, 1]:
               
                calcX = pixelSize * (imageCenter[0] + predictedPattern[predictedSpotRowIdx, 10] * numpy.cos(predictedPattern[predictedSpotRowIdx, 8])) ### Lab frame position, in m, relative to direct beam ###
                calcY = pixelSize * (imageCenter[1] + predictedPattern[predictedSpotRowIdx, 10] * numpy.sin(predictedPattern[predictedSpotRowIdx, 8])) 
                
                distanceFromPrediction = numpy.sqrt((calcX - detectedPeaks[detectedPeakRowIdx, 2])**2+(calcY - detectedPeaks[detectedPeakRowIdx, 3])**2)/pixelSize   
                detectedPeaks[detectedPeakRowIdx, 5] = distanceFromPrediction # Updated distance of detected peak from refined predicted peak, in pxls.
                break

    return detectedPeaks  