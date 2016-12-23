import numpy
cimport numpy
cimport cython

DTYPE = numpy.float32
ctypedef numpy.float32_t DTYPE_t

DTYPEI = numpy.int16
ctypedef numpy.int16_t DTYPEI_t

#@cython.boundscheck(False)
#@cython.wraparound(False)

def discontinuityMask(numpy.ndarray[dtype=DTYPE_t, ndim=2] xGeometry, numpy.ndarray[dtype=DTYPE_t, ndim=2] yGeometry,  
		      numpy.ndarray[dtype=numpy.int_t, ndim=2] modules):

    cdef int i
    cdef int j
    cdef float deltaRight_geoX
    cdef float deltaRight_geoY
    cdef float deltaTop_geoX
    cdef float deltaTop_geoY
    
    
    for i in range(0, xGeometry.shape[0]-1):
        for j in range(0, xGeometry.shape[1]-1):
            deltaRight_geoX = xGeometry[i, j+1] - xGeometry[i, j]
            deltaRight_geoY = yGeometry[i, j+1] - yGeometry[i, j]
            
            deltaTop_geoX = xGeometry[i+1, j] - xGeometry[i, j]
            deltaTop_geoY = yGeometry[i+1, j] - yGeometry[i, j]
            
            if abs(deltaRight_geoX) >= 0.000110 or abs(deltaRight_geoY) >= 0.000110 or abs(deltaTop_geoX) >= 0.000110 or abs(deltaTop_geoY) >= 0.000110:
                modules[i, j] = 1 
    
    return modules