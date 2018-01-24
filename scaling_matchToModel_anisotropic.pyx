cimport numpy
import numpy
cimport cython
import random
from libc.math cimport sqrt  # pow, sin, cos, fabs

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t


def matching(DTYPE_t[:, :] spotsL1, 
             DTYPE_t[:, :] spotsL2, 
             float deltaQrodThreshold, 
             float cellSize, 
             float resolution_limit_3D):
    # spotsL1: h k qRod I flag              L1 can be the model
    # spotsL2: h k qRod I flag

    cdef int n_pairs = 0
   
    cdef int nSpotsL1 = spotsL1.shape[0]
    cdef int nSpotsL2 = spotsL2.shape[0]
    
    cdef int index
    
    cdef float q_x
    cdef float q_y
    cdef float q_3D
    cdef float resolution_3D
        
    I1 = []
    I2 = []
    qs_2D = []
    qs_rod =[]
    
    ####################################
    #cdef float resolutionLimit = 7.1 # A
    #cdef float cellSize = 62.45      # A
    ####################################
    
    h_L1 = spotsL1[:, 0]
    k_L1 = spotsL1[:, 1]
    q_L1 = spotsL1[:, 2]
    I_L1 = spotsL1[:, 3]
    q_L1 = numpy.asarray(q_L1, dtype = numpy.float32)
    
    #################################################################################################
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    #################################################################################################
    
    for i in range(0, nSpotsL2):
        h = spotsL2[i, 0]
        k = spotsL2[i, 1]
        q = spotsL2[i, 2] # qRod
        I = spotsL2[i, 3]
        
        ############################################
        reciprocalVector = [h, k]*reciprocalCellRows
        q_x = reciprocalVector[0,0]                # A^(-1)
        q_y = reciprocalVector[0,1]                # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)         # A^(-1)
        q_3D = numpy.sqrt(q_x**2 + q_y**2 + q**2)  # A^(-1)
        resolution_3D = 2* numpy.pi / q_3D         # A
#        resolution = 2* numpy.pi / q_2D           # A
        ############################################
                
        if numpy.isnan(spotsL2[i, 3]):
            continue  
        
        ################################     
        #if resolution < resolutionLimit:
        if resolution_3D < resolution_limit_3D:
            continue
        ################################
        
        indices_qEqual    = numpy.argwhere(abs(q_L1 - q) <= deltaQrodThreshold)
        indices_qOpposite = numpy.argwhere(abs(q_L1 + q) <= deltaQrodThreshold)
        
        for index in indices_qEqual:
            if numpy.isnan(spotsL1[index, 3]):
                continue
            h1 = h_L1[index]
            k1 = k_L1[index]
            if ((h1 == h and k1 == k) or (h1 == -h-k and k1 == h) or (h1 == k and k1 == -h-k)):
                #if I_L1[index] > 0 and I > 0:
                I1.append(I_L1[index])
                I2.append(I)
                qs_2D.append(q_2D)
                qs_rod.append(q)
        for index in indices_qOpposite:
            if numpy.isnan(spotsL1[index, 3]):
                continue   
            h1 = h_L1[index]
            k1 = k_L1[index]
            if ((h1 == -h and k1 == -k) or (h1 == h+k and k1 == -h) or (h1 == -k and k1 == h+k)):
                #if I_L1[index] > 0 and I > 0:
                I1.append(I_L1[index])
                I2.append(I)
                qs_2D.append(q_2D)
                qs_rod.append(q)
       
    n_pairs = len(I1)                
    
    return n_pairs, I1, I2, qs_2D, qs_rod