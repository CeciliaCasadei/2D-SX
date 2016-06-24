import numpy
cimport numpy
cimport cython

DTYPE = numpy.float64
ctypedef numpy.float64_t DTYPE_t

DTYPEI = numpy.int16
ctypedef numpy.int16_t DTYPEI_t

@cython.boundscheck(False)
@cython.wraparound(False)

def unassembledMatching(numpy.ndarray[dtype=DTYPEI_t, ndim=2] xIndices, numpy.ndarray[dtype=DTYPEI_t, ndim=2] yIndices):
    cdef int nMatches = 0
    cdef int successFlag = 0
    cdef int myRowIndex_x
    cdef int myRowIndex_y
    cdef int i
    cdef int j
    cdef int i_good = 0
    cdef int j_good = 0
    cdef int myRowIndex_x_limit = xIndices.shape[0]
    cdef int myRowIndex_y_limit = yIndices.shape[0]
    
    for myRowIndex_x in range(0, myRowIndex_x_limit): 
        i = xIndices[myRowIndex_x, 0] # row index
        j = xIndices[myRowIndex_x, 1] # column index
        for myRowIndex_y in range(0, myRowIndex_y_limit): 
            if yIndices[myRowIndex_y, 0] == i and yIndices[myRowIndex_y, 1] == j:
                nMatches = nMatches + 1
                i_good = i
                j_good = j
                
    if nMatches == 1:  
        successFlag = 1
            
    return successFlag, i_good, j_good