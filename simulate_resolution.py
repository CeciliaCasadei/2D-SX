# -*- coding: utf-8 -*-
import numpy
cellSize = 62.45 #62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2* numpy.pi * directCell.I                                                   # A^(-1)

h = 5
k = 9
qRod = 0 #0.2991 # A^(-1)   
reciprocalVector = [h, k]*reciprocalCellRows
q_x = reciprocalVector[0,0]         # A^(-1)
q_y = reciprocalVector[0,1]         # A^(-1)
q = numpy.sqrt(q_x**2 + q_y**2 + qRod**2)     # A^(-1)
if q != 0:
    resolution = 2* numpy.pi / q    # A
print q_x
print q_y    
print resolution
print q