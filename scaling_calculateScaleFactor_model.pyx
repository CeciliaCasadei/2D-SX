cimport numpy
import numpy
cimport cython
import random
from libc.math cimport sqrt  # pow, sin, cos, fabs

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t


def calculateScaleFactorFunction(DTYPE_t[:, :] spotsModel, DTYPE_t[:, :] spotsL, float deltaQrodThreshold):
    # spotsModel: h k qRod I 
    # spotsL:     n h k qRod I Icorrected h_transformed k_transformed      

    cdef int n_pairs = 0
   
    cdef int nSpotsmodel = spotsModel.shape[0]
    cdef int nSpotsL     = spotsL.shape[0]
    
    cdef int index
    cdef DTYPE_t h, k  # why?
        
    I1 = []
    I2 = []
    
    hModel = spotsModel[:, 0]
    kModel = spotsModel[:, 1]
    qModel = spotsModel[:, 2]
    Imodel = spotsModel[:, 3]
    qModel = numpy.asarray(qModel, dtype = numpy.float32)
    
    for i in range(0, nSpotsL):
        h = spotsL[i, 6]
        k = spotsL[i, 7]
        q = spotsL[i, 3]
        I = spotsL[i, 5]
        
        if numpy.isnan(spotsL[i, 5]):
            continue     
        
        indices_qEqual    = numpy.argwhere(abs(qModel - q) <= deltaQrodThreshold)
        indices_qOpposite = numpy.argwhere(abs(qModel + q) <= deltaQrodThreshold)
        
        for index in indices_qEqual:
            if ((hModel[index] == h and kModel[index] == k) or (hModel[index] == -h-k and kModel[index] == h) or (hModel[index] == k and kModel[index] == -h-k)):
                I1.append(Imodel[index])
                I2.append(I)
        for index in indices_qOpposite:
            if ((hModel[index] == -h and kModel[index] == -k) or (hModel[index] == h+k and kModel[index] == -h) or (hModel[index] == -k and kModel[index] == h+k)):
                I1.append(Imodel[index])
                I2.append(I)
       
    n_pairs = len(I1)                
    if n_pairs > 5:
        I1 = numpy.asarray(I1)
        I2 = numpy.asarray(I2)
        I1 = I1[:,numpy.newaxis]

        scale, _, _, _ = numpy.linalg.lstsq(I1, I2) # I2 = scale*I1    Lattice = scale*Model   scale = scale_modelToLattice 
    else:
        scale = numpy.nan
    
    return n_pairs, scale

