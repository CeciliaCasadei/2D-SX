# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot
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
    
max_resolution = 7
    
Is_matlab = []
sigmas_matlab = []
IoverSigI_matlab = []

Is_python = []
sigmas_python = []
IoverSigI_python = []

cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I                                                  # A^(-1)

intensityFile_matlab = open('./Output_r0127/transformAndScale/MATLAB_OrbitIntensities_final_r0127_Indexable.txt', 'r')
intensityFile_python = open('./Output_r0127/transformAndScale/mergedIntensities_h_k_qRod_l_avgI_sigmaI.txt', 'r')
intensityFile_python_list = list(intensityFile_python)
nList = len(intensityFile_python_list)


### MATLAB - PYTHON ip TRANSFORMATION ###
n = 0
for I_line_matlab in intensityFile_matlab:
    h = int(I_line_matlab[2:5])
    k = int(I_line_matlab[6:9])
    
    reciprocalVector = [h, k]*reciprocalCellRows
    q_x = reciprocalVector[0,0]         # A^(-1)
    q_y = reciprocalVector[0,1]         # A^(-1)
    q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
    resolution = 2*numpy.pi/q
    if resolution < max_resolution:
        continue
   
    I = float(I_line_matlab[59:66])
    sigI = float(I_line_matlab[70:77])
    print '\n', h, k
    for i in range(0, nList):
        I_line_python = intensityFile_python_list[i]
        
        splittedLine = I_line_python.split()
        h_python = int(splittedLine[0])
        k_python = int(splittedLine[1])
        I_python = float(splittedLine[4])
        sigI_python = float(splittedLine[5])
        if k == 0:
            if h == h_python and k_python == 0:
                print 'matched', h, k, I, sigI, 'to', h_python, k_python, I_python, sigI_python
                Is_matlab.append(I)
                Is_python.append(I_python)
                sigmas_matlab.append(sigI)
                sigmas_python.append(sigI_python)
                IoverSigI_matlab.append(I/sigI)
                IoverSigI_python.append(I_python/sigI_python)
                n = n+1
        elif h == -k_python and k == -h_python:
            print 'matched', h, k, I, sigI, 'to', h_python, k_python, I_python, sigI_python
            Is_matlab.append(I)
            Is_python.append(I_python)
            sigmas_matlab.append(sigI)
            sigmas_python.append(sigI_python)
            IoverSigI_matlab.append(I/sigI)
            IoverSigI_python.append(I_python/sigI_python)
            n = n+1
            
print n

CC = Correlate(Is_matlab, Is_python)
print CC
print 'I no sum (circle) python - I no sum (circle) matlab (n pairs = %d): CC = %.4f'%(len(Is_python), CC)
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.4f'%(len(Is_python), CC))
matplotlib.pyplot.scatter(Is_python, Is_matlab)
matplotlib.pyplot.gca().set_xlabel('photon counts - no sum - integral on circle - python')
matplotlib.pyplot.gca().set_ylabel('photon counts - no sum - integral on circle - matlab', rotation = 'vertical')
matplotlib.pyplot.savefig('./Output_r0127/transformAndScale/Is_matlab_VS_Is_python_%dA.png'%max_resolution)
matplotlib.pyplot.close()

CC = Correlate(IoverSigI_matlab, IoverSigI_python)
print CC
print 'I/sigI, no sum (circle) python - I/sigI, no sum (circle) matlab (n pairs = %d): CC = %.4f'%(len(Is_python), CC)
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.4f'%(len(Is_python), CC))
matplotlib.pyplot.scatter(Is_python, Is_matlab)
matplotlib.pyplot.gca().set_xlabel('I/sigI - no sum - integral on circle - python')
matplotlib.pyplot.gca().set_ylabel('I/sigI - no sum - integral on circle - matlab', rotation = 'vertical')
matplotlib.pyplot.axes().set_aspect('equal', 'datalim')
matplotlib.pyplot.savefig('./Output_r0127/transformAndScale/IoverSigI_matlab_VS_IoverSigI_python_%dA.png'%max_resolution)
matplotlib.pyplot.close()
