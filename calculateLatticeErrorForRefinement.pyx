# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 18:05:55 2016

@author: casadei_c
"""
import numpy
cimport numpy
cimport cython

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)

def calculateLatticeError(numpy.ndarray[dtype=DTYPE_t, ndim=2] myPredictedPattern, numpy.ndarray[dtype=DTYPE_t, ndim=2] indexedPeaksTable):
                                                  
    cdef float latticeError = 0
    cdef int nMatchedPeaks = 0
    cdef float spotError
    cdef int indexedPeakRowIdx
    cdef int indexedPeaksNrows = indexedPeaksTable.shape[0]
    cdef int predictedSpotRowIdx
    cdef int predictedPatternNrows = myPredictedPattern.shape[0]
    cdef float obsX
    cdef float obsY
    cdef float calcX
    cdef float calcY
           
    for indexedPeakRowIdx in range(0, indexedPeaksNrows):
        obsX = indexedPeaksTable[indexedPeakRowIdx, 2] * numpy.cos(indexedPeaksTable[indexedPeakRowIdx, 3])
        obsY = indexedPeaksTable[indexedPeakRowIdx, 2] * numpy.sin(indexedPeaksTable[indexedPeakRowIdx, 3])
        for predictedSpotRowIdx in range(0, predictedPatternNrows):
            if myPredictedPattern[predictedSpotRowIdx, 0] == indexedPeaksTable[indexedPeakRowIdx, 0] and myPredictedPattern[predictedSpotRowIdx, 1] == indexedPeaksTable[indexedPeakRowIdx, 1]:
                nMatchedPeaks = nMatchedPeaks + 1
                calcX = myPredictedPattern[predictedSpotRowIdx, 10] * numpy.cos(myPredictedPattern[predictedSpotRowIdx, 8])
                calcY = myPredictedPattern[predictedSpotRowIdx, 10] * numpy.sin(myPredictedPattern[predictedSpotRowIdx, 8])
                spotError = numpy.sqrt((calcX-obsX)**2+(calcY-obsY)**2)/numpy.sqrt(obsX**2+obsY**2)
                latticeError = latticeError + spotError
                break
        
    return latticeError, nMatchedPeaks