cimport numpy
import numpy
cimport cython
import random
from libc.math cimport sqrt  # pow, sin, cos, fabs

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t


def calculateScaleFactorFunction(DTYPE_t[:, :] spotsL1, DTYPE_t[:, :] spotsL2, float deltaQrodThreshold):
    # spotsLi : n h k qRod I

    cdef int n_pairs = 0

    
    cdef int nSpotsL1 = spotsL1.shape[0]
    cdef int nSpotsL2 = spotsL2.shape[0]
    
    cdef int i, j
    cdef DTYPE_t h, k  # why?
        

    I1 = []
    I2 = []

    
    for i in range(0, nSpotsL1):
        h = spotsL1[i, 1]
        k = spotsL1[i, 2]
        
        if numpy.isnan(spotsL1[i, 4]):
            continue
        
        for j in range(0, nSpotsL2):
            if numpy.isnan(spotsL2[j, 4]):
                continue
            
            
            if (spotsL2[j, 1] == h and spotsL2[j, 2] == k) or (spotsL2[j, 1] == -h-k and spotsL2[j, 2] == h) or (spotsL2[j, 1] == k and spotsL2[j, 2] == -h-k):
                if abs(spotsL1[i, 3] - spotsL2[j, 3]) < deltaQrodThreshold:                
                    n_pairs = n_pairs + 1
                    I1.append(spotsL1[i, 4])
                    I2.append(spotsL2[j, 4])
            if (spotsL2[j, 1] == -h and spotsL2[j, 2] == -k) or (spotsL2[j, 1] == h+k and spotsL2[j, 2] == -h) or (spotsL2[j, 1] == -k and spotsL2[j, 2] == h+k):
                if abs(spotsL1[i, 3] + spotsL2[j, 3]) < deltaQrodThreshold:                
                    n_pairs = n_pairs + 1
                    I1.append(spotsL1[i, 4])
                    I2.append(spotsL2[j, 4])  
                    
    if n_pairs > 5:
        I1 = numpy.asarray(I1)
        I2 = numpy.asarray(I2)
        I1 = I1[:,numpy.newaxis]

        scale, _, _, _ = numpy.linalg.lstsq(I1, I2) # I2 = scale*I1  scale = scale_L1toL2
    else:
        scale = numpy.nan
    
    return n_pairs, scale

