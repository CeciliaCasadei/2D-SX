# -*- coding: utf-8 -*-
import numpy
def Correlate(x1, x2):
    x1Avg = numpy.average(x1)
    x2Avg = numpy.average(x2)
    numTerm = numpy.multiply(x1-x1Avg, x2-x2Avg)
    num = numTerm.sum()
    resX1Sq = numpy.multiply(x1-x1Avg, x1-x1Avg)
    resX2Sq = numpy.multiply(x2-x2Avg, x2-x2Avg)
    den = numpy.sqrt(numpy.multiply(resX1Sq.sum(), resX2Sq.sum()))
    CC = num/den
    return CC
    
    
Is_matlab = []
Is_python = []
intensityFile = open('./Output_r0127/transformAndScale/r0127_h_k_I_PYTHON_all_I_MATLAB_all.txt', 'r')
#intensityFile = open('./Output_r0127/transformAndScale/r0127_h_k_I_MATLAB_best10_I_PYTHON_all.txt', 'r')
for I_line in intensityFile:
    splittedLine = I_line.split()
    I_matlab = float(splittedLine[3])
    I_python = float(splittedLine[2])
    print '%f %f'%(I_matlab, I_python)
    Is_matlab.append(I_matlab)
    Is_python.append(I_python)
    
print 'MATLAB'
print Is_matlab
print "PYTHON"
print Is_python
print '\n'  
  
CC = Correlate(Is_matlab, Is_python)
print CC