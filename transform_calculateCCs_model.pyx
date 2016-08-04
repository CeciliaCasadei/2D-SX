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

        
    
def determineTransformation(DTYPE_t[:, :] spotsModel, DTYPE_t[:, :] spotsL, float deltaQrodThreshold):
    # spotsModel: h k qRod I
    # spotsLattice: n h k qRod I Icorrected

    cdef int n_I
    cdef int n_i
    cdef int n_p
    cdef int n_ip
    
    cdef int nSpotsModel = spotsModel.shape[0]
    cdef int nSpotsL = spotsL.shape[0]
    
    cdef int i, j
    cdef DTYPE_t h, k  # why?
        
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
    
    hModel = spotsModel[:, 0]
    kModel = spotsModel[:, 1]
    qModel = spotsModel[:, 2]
    Imodel = spotsModel[:, 3]
    
    for i in range(0, nSpotsL):
        h = spotsL[i, 1]
        k = spotsL[i, 2]
        q = spotsL[i, 3]
        I = spotsL[i, 5]
        
        if numpy.isnan(spotsL[i, 5]):
            continue
        
        qModel = numpy.asarray(qModel, dtype = numpy.float32)
        indices_qEqual    = numpy.argwhere(abs(qModel - q) <= deltaQrodThreshold)
        indices_qOpposite = numpy.argwhere(abs(qModel + q) <= deltaQrodThreshold)
        
        for index in indices_qEqual:
            # IDENTITY
            if ((hModel[index] == h and kModel[index] == k) or (hModel[index] == -h-k and kModel[index] == h) or (hModel[index] == k and kModel[index] == -h-k)):
                I1_I.append(Imodel[index])
                I2_I.append(I)
            # INVERSION
            if ((hModel[index] == -h and kModel[index] == -k) or (hModel[index] == h+k and kModel[index] == -h) or (hModel[index] == -k and kModel[index] == h+k)):
                I1_i.append(Imodel[index])
                I2_i.append(I)
            # PERMUTATION
            if ((hModel[index] == k and kModel[index] == h) or (hModel[index] == -k-h and kModel[index] == k) or (hModel[index] == h and kModel[index] == -k-h)):
                I1_p.append(Imodel[index])
                I2_p.append(I)
            # INVERSION - PERMUTATION
            if ((hModel[index] == -k and kModel[index] == -h) or (hModel[index] == k+h and kModel[index] == -k) or (hModel[index] == -h and kModel[index] == k+h)):
                I1_ip.append(Imodel[index])
                I2_ip.append(I)
                
        for index in indices_qOpposite:
            # INVERSION
            if ((hModel[index] == h and kModel[index] == k) or (hModel[index] == -h-k and kModel[index] == h) or (hModel[index] == k and kModel[index] == -h-k)):
                I1_i.append(Imodel[index])
                I2_i.append(I)
            # IDENTITY
            if ((hModel[index] == -h and kModel[index] == -k) or (hModel[index] == h+k and kModel[index] == -h) or (hModel[index] == -k and kModel[index] == h+k)):
                I1_I.append(Imodel[index])
                I2_I.append(I)
            # INVERSION-PERMUTATION
            if ((hModel[index] == k and kModel[index] == h) or (hModel[index] == -k-h and kModel[index] == k) or (hModel[index] == h and kModel[index] == -k-h)):
                I1_ip.append(Imodel[index])
                I2_ip.append(I)
            # PERMUTATION
            if ((hModel[index] == -k and kModel[index] == -h) or (hModel[index] == k+h and kModel[index] == -k) or (hModel[index] == -h and kModel[index] == k+h)):
                I1_p.append(Imodel[index])
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
    
    
    
    
##########################################################################################################
def DRAFT_determineTransformation(DTYPE_t[:, :] spotsModel, DTYPE_t[:, :] spotsL, float deltaQrodThreshold):
    # spotsModel: h k qRod I
    # spotsLattice: n h k qRod I Icorrected

    cdef int n_I
    cdef int n_i
    cdef int n_p
    cdef int n_ip
    
    cdef int nSpotsModel = spotsModel.shape[0]
    cdef int nSpotsL = spotsL.shape[0]
    
    cdef int i, j
    cdef DTYPE_t h, k  # why?
        
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
    
    hModel = spotsModel[:, 0]
    kModel = spotsModel[:, 1]
    qModel = spotsModel[:, 2]
    Imodel = spotsModel[:, 3]
    
    for i in range(0, nSpotsL):
        h = spotsL[i, 1]
        k = spotsL[i, 2]
        q = spotsL[i, 3]
        I = spotsL[i, 5]
        
        if numpy.isnan(spotsL[i, 5]):
            continue
        
        # IDENTITY
        Ipartial_model_I = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]-q)<=deltaQrodThreshold and (          (hModel[item] == h    and kModel[item] == k   ) 
                                                                                                                                  or (hModel[item] == -h-k and kModel[item] == h   )
                                                                                                                                  or (hModel[item] == k    and kModel[item] == -h-k)      )]
        for item in Ipartial_model_I:
            I1_I.append(item)
            I2_I.append(I)
        Ipartial_model_I = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]+q)<=deltaQrodThreshold and (          (hModel[item] == -h  and kModel[item] == -k   ) 
                                                                                                                                  or (hModel[item] == h+k and kModel[item] == -h   )  
                                                                                                                                  or (hModel[item] == -k  and kModel[item] == h+k  )      )]
        for item in Ipartial_model_I:
            I1_I.append(item)
            I2_I.append(I)
            
        # INVERSION
        Ipartial_model_i = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]-q)<=deltaQrodThreshold and (          (hModel[item] == -h   and kModel[item] == -k  ) 
                                                                                                                                  or (hModel[item] == h+k  and kModel[item] == -h  )
                                                                                                                                  or (hModel[item] == -k   and kModel[item] == h+k )      )]
        for item in Ipartial_model_i:
            I1_i.append(item)
            I2_i.append(I)
        Ipartial_model_i = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]+q)<=deltaQrodThreshold and (          (hModel[item] == h    and kModel[item] == k   ) 
                                                                                                                                  or (hModel[item] == -h-k and kModel[item] == h   ) 
                                                                                                                                  or (hModel[item] == k    and kModel[item] == -h-k)      )]
        for item in Ipartial_model_i:
            I1_i.append(item)
            I2_i.append(I)
            
        if h != k:  
            # PERMUTATION
            Ipartial_model_p = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]-q)<=deltaQrodThreshold  and (     (hModel[item] == k    and kModel[item] == h   ) 
                                                                                                                                  or (hModel[item] == -k-h and kModel[item] == k   ) 
                                                                                                                                  or (hModel[item] == h    and kModel[item] == -k-h)      )]
            for item in Ipartial_model_p:
                I1_p.append(item)
                I2_p.append(I)
            Ipartial_model_p = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]+q)<=deltaQrodThreshold  and (     (hModel[item] == -k   and kModel[item] == -h  ) 
                                                                                                                                  or (hModel[item] == k+h  and kModel[item] == -k  ) 
                                                                                                                                  or (hModel[item] == -h   and kModel[item] == k+h )      )]
            for item in Ipartial_model_p:
                I1_p.append(item)
                I2_p.append(I)
            # INVERSION-PERMUTATION
            Ipartial_model_ip = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]-q)<=deltaQrodThreshold and (     (hModel[item] == -k   and kModel[item] == -h  ) 
                                                                                                                                  or (hModel[item] == k+h  and kModel[item] == -k  ) 
                                                                                                                                  or (hModel[item] == -h   and kModel[item] == k+h )      )]
            for item in Ipartial_model_ip:
                I1_ip.append(item)
                I2_ip.append(I)
            Ipartial_model_ip = [Imodel[item] for item in range(0, len(Imodel)) if abs(qModel[item]+q)<=deltaQrodThreshold and (     (hModel[item] == k    and kModel[item] == h   ) 
                                                                                                                                  or (hModel[item] == -k-h and kModel[item] == k   ) 
                                                                                                                                  or (hModel[item] == h    and kModel[item] == -k-h)      )]
            for item in Ipartial_model_ip:
                I1_ip.append(item)
                I2_ip.append(I)

                        
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