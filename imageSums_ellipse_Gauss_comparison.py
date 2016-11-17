# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import scipy.optimize

import imageSums_utilities

outputFolder = './Output_imageSums_sigmaFits'

# RECIPROCAL CELL
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I             

# COMPARE I_SUM_ELLIPSE_PY, I_GAUSS_FIXED_SIGMAS_PY in RESOLUTION RANGE 55 - 4 A (AFTER FILTERING OUT BAD GAUSS FITS)
Is_gauss_fixed_sigmas_python = []
Is_sum_ellipse_python = []
qs = []

intensityFile_sum = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse.txt'%outputFolder, 'r')                                  # 55 - 4 A
for I_line in intensityFile_sum:   # 55 - 4 A
    splittedLine = I_line.split()
    h = int(splittedLine[0])
    k = int(splittedLine[1])
    I_gauss_fixed_sigmas_python = float(splittedLine[3])
    I_sum_ellipse_python        = float(splittedLine[4])
    
    reciprocalVector = [h, k]*reciprocalCellRows
    q_x = reciprocalVector[0,0]         # A^(-1)
    q_y = reciprocalVector[0,1]         # A^(-1)
    q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
    
    if not numpy.isnan(I_gauss_fixed_sigmas_python):
        Is_gauss_fixed_sigmas_python.append(I_gauss_fixed_sigmas_python)
        Is_sum_ellipse_python.append(I_sum_ellipse_python)
        qs.append(q)
    else:
        print 'I_gauss_fixed_sigmas_python isnan:', h, k
          
intensityFile_sum.close()

# LINE FIT
popt_sigY_line, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
offset = popt_sigY_line[0]
slope  = popt_sigY_line[1]
print 'Line fit: I_gauss = %.2f + %.2f * I_ellipse'%(offset, slope)

# CC CALCULATION
CC = imageSums_utilities.Correlate(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
print 'I sum (ellipse) python - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_sum_ellipse_python), CC)

# PLOTS
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python), CC), y=1.05)
matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
matplotlib.pyplot.gca().set_xscale('log')
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_xlim([0.01,1000])
matplotlib.pyplot.gca().set_ylim([0.01,1000])
matplotlib.pyplot.gca().set_xlabel(r'I$_{\rm sum, ellipse}$ (photon counts)')
matplotlib.pyplot.gca().set_ylabel(r'I$_{\rm sum, gaussian}$ (photon counts)')
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig('%s/paper_Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_logscale.png'%outputFolder, dpi=2*96)
matplotlib.pyplot.close()

matplotlib.pyplot.figure()
matplotlib.pyplot.scatter(qs, Is_gauss_fixed_sigmas_python, marker='o', edgecolor='none', c='b', s=8, label='Gaussian')
matplotlib.pyplot.scatter(qs, Is_sum_ellipse_python, c='m', marker='x', linewidth='0.4', s=10, label='Ellipse')
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_ylim([0.01,1000])
matplotlib.pyplot.gca().set_xlabel(r'q ($\AA^{-1}$)')
matplotlib.pyplot.gca().set_ylabel(r'I (photon counts)')
legend = matplotlib.pyplot.gca().legend(scatterpoints=1, fontsize=10)
legend.draw_frame(False)
matplotlib.pyplot.savefig('%s/paper_Is_gauss_fixed_sigmas_python_Is_sum_ellipse_python_Vs_q_4A.png'%outputFolder, dpi=2*96)
matplotlib.pyplot.close()