# cython: profile=False
# cython: linetrace=False
# cython: binding=False
# cython: boundscheck=False
# cython: wraparound=False

cimport numpy
import numpy
cimport cython
import random
from libc.math cimport sqrt  # pow, sin, cos, fabs

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t




cdef average(DTYPE_t[:] x):
    cdef int i
    cdef DTYPE_t total_sum    
    total_sum = 0

    if x.shape[0] == 0:
        print "[average]: cannot average on a zero-elements array, returning 0"
        return 0

    for i in range(x.shape[0]):
        total_sum += x[i]
    
    return total_sum / float(x.shape[0])



cdef sum_float(DTYPE_t[:] x):
    cdef int i
    cdef DTYPE_t total_sum
    
    total_sum = 0
        
    for i in range(x.shape[0]):
        total_sum += x[i]
    return total_sum
    
    
    
cdef subtract_multiply_float(DTYPE_t[:] x1, DTYPE_t x1_0, DTYPE_t[:] x2, DTYPE_t x2_0):
    cdef int i
    cdef DTYPE_t[:] total = numpy.zeros((len(x1)), dtype=numpy.float32)
    
    for i in range(x1.shape[0]):
        total[i] = (x1[i] - x1_0) * (x2[i] - x2_0)
    
    return total



    
cdef Correlate(DTYPE_t[:] x1, DTYPE_t[:] x2):
    cdef DTYPE_t x1Avg, x2Avg, num, den
    cdef DTYPE_t[:] numTerm, resX1Sq, resX2Sq
    
    x1Avg = average(x1)
    x2Avg = average(x2)
    
    numTerm = subtract_multiply_float(x1, x1Avg, x2, x2Avg)
    num = sum_float(numTerm)
    resX1Sq = subtract_multiply_float(x1, x1Avg, x1, x1Avg)
    resX2Sq = subtract_multiply_float(x2, x2Avg, x2, x2Avg)

    den = sqrt(sum_float(resX1Sq) * sum_float(resX2Sq))
    try:
        CC = num / den

    except:
        print "-----------------------------------------------------"
        print "ERROR in Correlate, cannot divide ", num, "by", den
        print "resX1Sq, resX1Sq.sum(), resX2Sq, resX2Sq.sum():", resX1Sq[:], sum_float(resX1Sq), resX2Sq[:], sum_float(resX2Sq)
        print "Returning CC=0"
        print "-----------------------------------------------------"
        return 0

    return CC



cdef sample_and_correlate(int n, int nmin, I1, I2, int size=100):
    cdef int k = 0
    cdef DTYPE_t[:] I1_extract, I2_extract, CC
    cdef DTYPE_t CC_avg, CC_item

    I1_extract = numpy.zeros((nmin), dtype=numpy.float32)
    I2_extract = numpy.zeros((nmin), dtype=numpy.float32)
    CC = numpy.zeros((size), dtype=numpy.float32)

    for i in range(size):
        k = 0
        mySample = random.sample(range(n), nmin)

        for j in mySample:
            I1_extract[k] = I1[j]
            I2_extract[k] = I2[j]
            k += 1
            
        CC_item = Correlate(I1_extract, I2_extract)
        CC[i] = CC_item
        
    CC_avg = average(CC)
    return CC_avg

        
    
def determineTransformation(DTYPE_t[:, :] spotsL1, DTYPE_t[:, :] spotsL2, float deltaQrodThreshold):
    # spotsL1: h k qRod I (flag)          # L1 CAN BE THE MODEL
    # spotsL2: h k qRod I flag

    cdef int n_I
    cdef int n_i
    cdef int n_p
    cdef int n_ip
    
    cdef int nSpotsL1 = spotsL1.shape[0]  # L1 or MODEL
    cdef int nSpotsL2 = spotsL2.shape[0]
    
    cdef int i, j
    cdef DTYPE_t h, k  # why?
    
    ####################################
    cdef float resolutionLimit = 7.1 # A
    cdef float cellSize = 62.45      # A
    ####################################
        
    n_I = 0
    I1_I = []
    I2_I = []
    n_i = 0
    I1_i = []
    I2_i = []
    n_p = 0
    I1_p = []
    I2_p = []
    n_ip = 0
    I1_ip = []
    I2_ip = []
    
    #################################################################################################
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    #################################################################################################
    
    h_L1 = spotsL1[:, 0]                  # L1 or MODEL
    k_L1 = spotsL1[:, 1]
    q_L1 = spotsL1[:, 2]
    I_L1 = spotsL1[:, 3]
    q_L1 = numpy.asarray(q_L1, dtype = numpy.float32)
    
    for i in range(0, nSpotsL2):
        h = spotsL2[i, 0]
        k = spotsL2[i, 1]
        q = spotsL2[i, 2]  # qRod
        I = spotsL2[i, 3]
        
        ############################################
        reciprocalVector = [h, k]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        resolution = 2* numpy.pi / q_2D     # A
        ############################################
        
        if numpy.isnan(spotsL2[i, 3]):
            continue
                
        ################################     
        if resolution < resolutionLimit:
            continue
        ################################
                
        indices_qEqual    = numpy.argwhere(abs(q_L1 - q) <= deltaQrodThreshold)
        indices_qOpposite = numpy.argwhere(abs(q_L1 + q) <= deltaQrodThreshold)
        
        for index in indices_qEqual:
            if numpy.isnan(spotsL1[index, 3]):
                continue
            h1 = h_L1[index]
            k1 = k_L1[index]
            # IDENTITY
            if ((h1 == h and k1 == k) or (h1 == -h-k and k1 == h) or (h1 == k and k1 == -h-k)):
                I1_I.append(I_L1[index])
                I2_I.append(I)
            # INVERSION
            if ((h1 == -h and k1 == -k) or (h1 == h+k and k1 == -h) or (h1 == -k and k1 == h+k)):
                I1_i.append(I_L1[index])
                I2_i.append(I)
            # PERMUTATION
            if ((h1 == k and k1 == h) or (h1 == -k-h and k1 == k) or (h1 == h and k1 == -k-h)):
                I1_p.append(I_L1[index])
                I2_p.append(I)
            # INVERSION - PERMUTATION
            if ((h1 == -k and k1 == -h) or (h1 == k+h and k1 == -k) or (h1 == -h and k1 == k+h)):
                I1_ip.append(I_L1[index])
                I2_ip.append(I)
                
        for index in indices_qOpposite:
            if numpy.isnan(spotsL1[index, 3]):
                continue
            h1 = h_L1[index]
            k1 = k_L1[index]
            # INVERSION
            if ((h1 == h and k1 == k) or (h1 == -h-k and k1 == h) or (h1 == k and k1 == -h-k)):
                I1_i.append(I_L1[index])
                I2_i.append(I)
            # IDENTITY
            if ((h1 == -h and k1 == -k) or (h1 == h+k and k1 == -h) or (h1 == -k and k1 == h+k)):
                I1_I.append(I_L1[index])
                I2_I.append(I)
            # INVERSION-PERMUTATION
            if ((h1 == k and k1 == h) or (h1 == -k-h and k1 == k) or (h1 == h and k1 == -k-h)):
                I1_ip.append(I_L1[index])
                I2_ip.append(I)
            # PERMUTATION
            if ((h1 == -k and k1 == -h) or (h1 == k+h and k1 == -k) or (h1 == -h and k1 == k+h)):
                I1_p.append(I_L1[index])
                I2_p.append(I)

                        
    n_I = len(I1_I)    
    n_i = len(I1_i)
    n_p = len(I1_p)
    n_ip = len(I1_ip)                                          
    n_min = min([n_I, n_i, n_p, n_ip])
    if n_min > 3:
        CC_I_avg = sample_and_correlate(n_I, n_min, I1_I, I2_I)
        CC_i_avg = sample_and_correlate(n_i, n_min, I1_i, I2_i)
        CC_p_avg = sample_and_correlate(n_p, n_min, I1_p, I2_p)
        CC_ip_avg = sample_and_correlate(n_ip, n_min, I1_ip, I2_ip)    
        
        avg_CCs = [CC_I_avg, CC_i_avg, CC_p_avg, CC_ip_avg]
    else:
        avg_CCs = [0, 0, 0, 0]
    return n_min, avg_CCs