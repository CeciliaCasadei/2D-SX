# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 10:09:03 2016
@author: casadei_c
CALCULATE RECIPROCAL LATTICE IN REFERENCE ORIENTATION:
[h k qx qy resolution(A) q]
IN ORDER OF ASCENDING q.
"""
import numpy

def buildReciprocalLatticeFunction(cellSize, hmax, kmax, resolutionLimit):
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I                                                   # A^(-1)
    reciprocalLattice = []
    i = 0
    for h in range(-hmax, +hmax+1):
        for k in range(-kmax, +kmax+1):
            reciprocalVector = [h, k]*reciprocalCellRows
            q_x = reciprocalVector[0,0]         # A^(-1)
            q_y = reciprocalVector[0,1]         # A^(-1)
            q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
            if q != 0:
                resolution = 2* numpy.pi / q    # A
            if q != 0 and resolution >= resolutionLimit:
                reciprocalLatticeItem = []
                reciprocalLatticeItem.append(int(h)) 
                reciprocalLatticeItem.append(int(k))
                reciprocalLatticeItem.append(q_x)
                reciprocalLatticeItem.append(q_y)
                reciprocalLatticeItem.append(resolution)
                reciprocalLatticeItem.append(q)
                reciprocalLattice.append(reciprocalLatticeItem)
                i = i + 1
    
    # SORTING REFLECTIONS IN ORDER OF ASCENDING q:
    for i_sort in range(0, i):
        for j_sort in range(i_sort+1, i):
            if reciprocalLattice[i_sort][5] > reciprocalLattice[j_sort][5]:
                h = reciprocalLattice[i_sort][0]
                k = reciprocalLattice[i_sort][1]
                qx = reciprocalLattice[i_sort][2]
                qy = reciprocalLattice[i_sort][3]
                resolution = reciprocalLattice[i_sort][4]
                q = reciprocalLattice[i_sort][5]
                reciprocalLattice[i_sort][0] = reciprocalLattice[j_sort][0]
                reciprocalLattice[i_sort][1] = reciprocalLattice[j_sort][1]
                reciprocalLattice[i_sort][2] = reciprocalLattice[j_sort][2]
                reciprocalLattice[i_sort][3] = reciprocalLattice[j_sort][3]
                reciprocalLattice[i_sort][4] = reciprocalLattice[j_sort][4]
                reciprocalLattice[i_sort][5] = reciprocalLattice[j_sort][5]
                reciprocalLattice[j_sort][0] = h
                reciprocalLattice[j_sort][1] = k
                reciprocalLattice[j_sort][2] = qx
                reciprocalLattice[j_sort][3] = qy
                reciprocalLattice[j_sort][4] = resolution
                reciprocalLattice[j_sort][5] = q
    # END SORTING
    
    return reciprocalLattice