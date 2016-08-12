cimport numpy
import numpy
cimport cython
import random
from libc.math cimport sqrt  # pow, sin, cos, fabs
import matplotlib.pyplot

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t


def calculateScaleFactorFunction(DTYPE_t[:, :] spotsL1, DTYPE_t[:, :] spotsL2, float deltaQrodThreshold):
    # spotsL1: h k qRod I flag              L1 can be the model
    # spotsL2: h k qRod I flag

    cdef int n_pairs = 0
   
    cdef int nSpotsL1 = spotsL1.shape[0]
    cdef int nSpotsL2 = spotsL2.shape[0]
    
    cdef int index
        
    I1 = []
    I2 = []
    
    h_L1 = spotsL1[:, 0]
    k_L1 = spotsL1[:, 1]
    q_L1 = spotsL1[:, 2]
    I_L1 = spotsL1[:, 3]
    q_L1 = numpy.asarray(q_L1, dtype = numpy.float32)
    
    for i in range(0, nSpotsL2):
        h = spotsL2[i, 0]
        k = spotsL2[i, 1]
        q = spotsL2[i, 2]
        I = spotsL2[i, 3]
        
        if numpy.isnan(spotsL2[i, 3]):
            continue     
        
        indices_qEqual    = numpy.argwhere(abs(q_L1 - q) <= deltaQrodThreshold)
        indices_qOpposite = numpy.argwhere(abs(q_L1 + q) <= deltaQrodThreshold)
        
        for index in indices_qEqual:
            if numpy.isnan(spotsL1[index, 3]):
                continue
            h1 = h_L1[index]
            k1 = k_L1[index]
            if ((h1 == h and k1 == k) or (h1 == -h-k and k1 == h) or (h1 == k and k1 == -h-k)):
                I1.append(I_L1[index])
                I2.append(I)
        for index in indices_qOpposite:
            if numpy.isnan(spotsL1[index, 3]):
                continue   
            h1 = h_L1[index]
            k1 = k_L1[index]
            if ((h1 == -h and k1 == -k) or (h1 == h+k and k1 == -h) or (h1 == -k and k1 == h+k)):
                I1.append(I_L1[index])
                I2.append(I)
       
    n_pairs = len(I1)                
    if n_pairs > 5:
        I1 = numpy.asarray(I1)
        I2 = numpy.asarray(I2)
        I1_toFit = I1[:,numpy.newaxis]

        scale, _, _, _ = numpy.linalg.lstsq(I1_toFit, I2) # I2 = scale*I1    Lattice2 = scale*Lattice1   scale = scale_L1toL2
        
    else:
        scale = numpy.nan
    
    return n_pairs, scale, I1, I2