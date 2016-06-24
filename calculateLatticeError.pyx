import numpy
cimport numpy
cimport cython

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def calculateMatrixElement(numpy.ndarray[dtype=DTYPE_t, ndim=1] imageCenter,
                           numpy.ndarray[dtype=DTYPE_t, ndim=2] predictedPattern, 
                           numpy.ndarray[dtype=DTYPE_t, ndim=2] detectedPeaks,
                           float pixelSize, int distanceThreshold):
                            
    cdef float latticeError = 0
    cdef float spotError
    cdef int detectedPeakRowIdx
    cdef int detectedPeaksNrows = detectedPeaks.shape[0]
    cdef int predictedSpotRowIdx
    cdef int predictedPatternNrows = predictedPattern.shape[0]
    cdef float calcX
    cdef float calcY
    cdef int nMatchedPeaks = 0
       
    for detectedPeakRowIdx in range(0, detectedPeaksNrows):
        if detectedPeaks[detectedPeakRowIdx, 5] <= distanceThreshold:
            for predictedSpotRowIdx in range(0, predictedPatternNrows):
                if predictedPattern[predictedSpotRowIdx, 0] == detectedPeaks[detectedPeakRowIdx, 0] and predictedPattern[predictedSpotRowIdx, 1] == detectedPeaks[detectedPeakRowIdx, 1]:
                    nMatchedPeaks = nMatchedPeaks + 1
                    calcX = pixelSize * (imageCenter[0] + predictedPattern[predictedSpotRowIdx, 10] * numpy.cos(predictedPattern[predictedSpotRowIdx, 8])) ### Lab frame position, in m, relative to direct beam ###
                    calcY = pixelSize * (imageCenter[1] + predictedPattern[predictedSpotRowIdx, 10] * numpy.sin(predictedPattern[predictedSpotRowIdx, 8])) 
                    
                    spotError = numpy.sqrt((calcX - detectedPeaks[detectedPeakRowIdx, 2])**2 + (calcY - detectedPeaks[detectedPeakRowIdx, 3])**2)/numpy.sqrt(calcX**2 + calcY**2)                    
                    latticeError = latticeError + spotError
                    break

    return nMatchedPeaks, latticeError