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
		      numpy.ndarray[dtype=numpy.int_t, ndim=2] maskBoxMatrix, int i_good, int j_good, int boxWidth):

    cdef int increment = 0
    cdef int leftBorder
    cdef int rightBorder
    cdef int topBorder
    cdef int bottomBorder
    cdef int leftBorderTruncation = -1
    cdef int rightBorderTruncation = -1
    cdef int topBorderTruncation = -1
    cdef int bottomBorderTruncation = -1
    cdef float deltaLeft_geoX
    cdef float deltaLeft_geoY
    cdef float deltaRight_geoX
    cdef float deltaRight_geoY
    cdef float deltaBottom_geoX
    cdef float deltaBottom_geoY
    cdef float deltaTop_geoX
    cdef float deltaTop_geoY
    cdef int maskBorder_left = 0
    cdef int maskBorder_down = 0

    ### SCAN ROW TO THE LEFT ###
    leftBorder = boxWidth
    for increment in range(0, boxWidth):
        if j_good - (increment+1) >= 0:
            deltaLeft_geoX = xGeometry[i_good, j_good-increment] - xGeometry[i_good, j_good-(increment+1)]
            deltaLeft_geoY = yGeometry[i_good, j_good-increment] - yGeometry[i_good, j_good-(increment+1)]
            if abs(deltaLeft_geoX) >= 0.000110 or abs(deltaLeft_geoY) >= 0.000110:
                leftBorder = increment
                maskBoxMatrix[:,0:(boxWidth-leftBorder+1)] = 0 # Columns from 0 to boxWidth-leftBorder included are set to zero
                maskBorder_left = boxWidth-leftBorder+1
                break
        else:
            leftBorderTruncation = increment
            break
    ### SCAN ROW TO THE RIGHT ###
    rightBorder = boxWidth
    for increment in range(0, boxWidth):
        if j_good + (increment+1) < xGeometry.shape[1]:
            deltaRight_geoX = xGeometry[i_good, j_good+increment] - xGeometry[i_good, j_good+(increment+1)]
            deltaRight_geoY = yGeometry[i_good, j_good+increment] - yGeometry[i_good, j_good+(increment+1)]
            if abs(deltaRight_geoX) >= 0.000110 or abs(deltaRight_geoY) >= 0.000110:
                rightBorder = increment
                maskBoxMatrix[:,(boxWidth+rightBorder):] = 0
                break
        else:
            rightBorderTruncation = increment
            break
    ### SCAN COLUMN TO THE BOTTOM ###
    bottomBorder = boxWidth
    for increment in range(0, boxWidth):
        if i_good - (increment+1) >= 0:
            deltaBottom_geoX = xGeometry[i_good-increment, j_good] - xGeometry[i_good-(increment+1), j_good]
            deltaBottom_geoY = yGeometry[i_good-increment, j_good] - yGeometry[i_good-(increment+1), j_good]
            if abs(deltaBottom_geoX) >= 0.000110 or abs(deltaBottom_geoY) >= 0.000110:
                bottomBorder = increment
                maskBoxMatrix[0:(boxWidth-bottomBorder+1),:] = 0 # Rows from 0 to boxWidth-bottomBorder included are set to zero
                maskBorder_down = boxWidth-bottomBorder+1
                break
        else:
            bottomBorderTruncation = increment
            break
    ### SCAN COLUMN TO THE TOP ###
    topBorder = boxWidth
    for increment in range(0, boxWidth):
        if i_good + (increment+1) < xGeometry.shape[0]:
            deltaTop_geoX = xGeometry[i_good+increment, j_good] - xGeometry[i_good+(increment+1), j_good]
            deltaTop_geoY = yGeometry[i_good+increment, j_good] - yGeometry[i_good+(increment+1), j_good]
            if abs(deltaTop_geoX) >= 0.000110 or abs(deltaTop_geoY) >= 0.000110:
                topBorder = increment
                maskBoxMatrix[(boxWidth+topBorder):,:] = 0
                break
        else:
            topBorderTruncation = increment
            break

    if leftBorderTruncation != -1:
        #print 'Left border truncation'
        maskBoxMatrix = maskBoxMatrix[:,(boxWidth-leftBorderTruncation):]
    if rightBorderTruncation != -1:
        #print 'Right border truncation'
        maskBoxMatrix =  maskBoxMatrix[:,0:(boxWidth+rightBorderTruncation+1)]
    if bottomBorderTruncation != -1:
        #print 'Bottom border truncation'
        maskBoxMatrix =  maskBoxMatrix[(boxWidth-bottomBorderTruncation):,:]
    if topBorderTruncation != -1:
        #print 'Top border truncation'
        maskBoxMatrix =  maskBoxMatrix[0:(boxWidth+topBorderTruncation+1),:]
    
    return maskBorder_left, maskBorder_down, maskBoxMatrix   
