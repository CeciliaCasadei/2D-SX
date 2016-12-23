# -*- coding: utf-8 -*-
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import scipy.optimize
import pickle

import imageSums_utilities


# PARAMETERS
photonThreshold = 2
filter_distanceThreshold = 2.0
nCountsPerPhoton = 26
ellipse_multiplicative_factor = 2.5 

# TO BE FITTED
qs = []
sigXs= []
sigYs= []

Is = []

# FOLDERS
outputFolder = './Output_imageSums_sigmaFits_compare_fits'
if not os.path.exists('%s'%outputFolder):
    os.mkdir('%s'%outputFolder)

# RECIPROCAL CELL
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I                                                  # A^(-1)

# LOOP ON SPOTS WITH GAUSS INTENSITY ABOVE 1 PHOTON AND GOOD GAUSSIAN FIT
fOpen = open('./Output_imageSums/h_k_Isum_Igauss_sigX_sigY.txt', 'r')
for hk_line in fOpen:
    splittedLine = hk_line.split()
    h      = int(splittedLine[0])
    k      = int(splittedLine[1])
    Isum   = float(splittedLine[2])
    Igauss = float(splittedLine[3])
    sigX   = float(splittedLine[4])
    sigY   = float(splittedLine[5])
    
    reciprocalVector = [h, k]*reciprocalCellRows
    q_x = reciprocalVector[0,0]         # A^(-1)
    q_y = reciprocalVector[0,1]         # A^(-1)
    q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
    
    if Igauss < photonThreshold*nCountsPerPhoton:
        print '(%4d, %4d) DISCARDED: I < %f photons'%(h, k, photonThreshold)
        continue
    if sigX > 7 or sigY > 7:
        print '(%4d, %4d) DISCARDED: Gauss width is too large.'%(h, k)
        continue

    qs.append(q)
    sigXs.append(sigX)
    sigYs.append(sigY)
    Is.append(Isum)

### FIND SIGMA VS Q BEHAVIOUR ###

# FIT sigX VS q WITH A QUADRATIC
popt_sigX_quadratic, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic, qs, sigXs)
x = numpy.linspace(0.1, 1.6, 100)
y_quadratic = imageSums_utilities.quadratic(x, *popt_sigX_quadratic)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigX VS q WITH A QUADRATIC WITH NO FIRST ORDER TERM
popt_sigX_quadratic_no_first_order, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic_no_first_order, qs, sigXs)
y_quadratic_no_first_order = imageSums_utilities.quadratic_no_first_order(x, *popt_sigX_quadratic_no_first_order)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_no_first_order, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic_no_first_order)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit_no_first_order.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigX VS q WITH A QUADRATIC AT LOW q
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.7]
low_q_sigXs = [sigXs[i] for i in range(0, len(qs)) if qs[i] < 0.7]
popt_sigX_quadratic_low_q, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic, low_qs, low_q_sigXs)
y_quadratic_low_q = imageSums_utilities.quadratic(x, *popt_sigX_quadratic_low_q)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_low_q, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic_low_q)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit_low_q.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigX VS q WITH A QUADRATIC AT LOW q, EXCLUDE FIRST POINT
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.7 and qs[i] > 0.15]
low_q_sigXs = [sigXs[i] for i in range(0, len(qs)) if qs[i] < 0.7 and qs[i] > 0.15]
popt_sigX_quadratic_low_q, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic, low_qs, low_q_sigXs)
y_quadratic_low_q = imageSums_utilities.quadratic(x, *popt_sigX_quadratic_low_q)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_low_q, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic_low_q)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit_low_q_exclude_1st_point.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigX VS q WITH A QUADRATIC WITH NO FIRST ORDER TERM, AT LOW q, EXCLUDE FIRST POINT
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.7 and qs[i] > 0.15]
low_q_sigXs = [sigXs[i] for i in range(0, len(qs)) if qs[i] < 0.7 and qs[i] > 0.15]
popt_sigX_quadratic_no_first_order_low_q, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic_no_first_order, low_qs, low_q_sigXs)
y_quadratic_no_first_order_low_q = imageSums_utilities.quadratic_no_first_order(x, *popt_sigX_quadratic_no_first_order_low_q)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_no_first_order_low_q, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic_no_first_order_low_q)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit_no_first_order_low_q_exclude_1st_point.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigX VS q WITH A QUADRATIC WITH NO FIRST ORDER TERM, AT LOW q
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.7]
low_q_sigXs = [sigXs[i] for i in range(0, len(qs)) if qs[i] < 0.7]
popt_sigX_quadratic_no_first_order_low_q, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic_no_first_order, low_qs, low_q_sigXs)
y_quadratic_no_first_order_low_q = imageSums_utilities.quadratic_no_first_order(x, *popt_sigX_quadratic_no_first_order_low_q)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_no_first_order_low_q, label='fit')
matplotlib.pyplot.title('%s'%popt_sigX_quadratic_no_first_order_low_q)
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit_no_first_order_low_q.png'%outputFolder)
matplotlib.pyplot.close()

offset = popt_sigX_quadratic_no_first_order_low_q[0]
secondOrder = popt_sigX_quadratic_no_first_order_low_q[1]

# FIT sigX VS q WITH A QUADRATIC WITH NO FIRST ORDER + SIGMA:
popt_sigX, pcov = scipy.optimize.curve_fit(lambda x, x0, k, scale: imageSums_utilities.quadratic_plus_sigmoid(x, offset, secondOrder, x0, k, scale), qs, sigXs)
y_quadratic_plus_sigmoid = imageSums_utilities.quadratic_plus_sigmoid(x, offset, secondOrder, popt_sigX[0], popt_sigX[1], popt_sigX[2])

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic_plus_sigmoid, label='fit')
matplotlib.pyplot.plot(x, y_quadratic_no_first_order_low_q)
matplotlib.pyplot.title('%s, %s'%(popt_sigX_quadratic_no_first_order_low_q, popt_sigX))
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticPlusSigmoidFit.png'%outputFolder)
matplotlib.pyplot.close()

# FIT sigY VS q WITH A LINE AT LOW q
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.9]
low_q_sigYs = [sigYs[i] for i in range(0, len(qs)) if qs[i] < 0.9]

popt_sigY_line, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, low_qs, low_q_sigYs)
y_line = imageSums_utilities.line(x, *popt_sigY_line)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_line, label='fit')
matplotlib.pyplot.title('%s'%popt_sigY_line)
matplotlib.pyplot.savefig('%s/sigYvsQ_lineFit.png'%outputFolder)
matplotlib.pyplot.close()

offset = popt_sigY_line[0]
slope  = popt_sigY_line[1]

# FIT sigY VS q WITH A LINE + SIGMA
popt_sigY, pcov = scipy.optimize.curve_fit(lambda x, x0, k, scale: imageSums_utilities.line_plus_sigmoid(x, offset, slope, x0, k, scale), qs, sigYs)
y_line_plus_sigmoid = imageSums_utilities.line_plus_sigmoid(x, offset, slope, popt_sigY[0], popt_sigY[1], popt_sigY[2])

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_line_plus_sigmoid, label='fit')
matplotlib.pyplot.plot(x, y_line)
matplotlib.pyplot.title('%s, %s'%(popt_sigY_line, popt_sigY))
matplotlib.pyplot.savefig('%s/sigYvsQ_linePlusSigmidFit.png'%outputFolder)
matplotlib.pyplot.close()

### SUMMARY PLOT ###
myFigure, (ax1, ax2) = matplotlib.pyplot.subplots(2, 1)    
ax1.scatter(qs, sigXs  , color='b', label='data', alpha=0.5)
ax1.set_title('Radial width', fontsize=6)
ax1.set_ylabel(r'$\sigma_{\rm radial}$ (pixels)')
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.plot(x, y_quadratic_no_first_order_low_q, label='fit', color='m', linestyle=':')
ax1.plot(x, y_quadratic_plus_sigmoid, label='fit', color='m')

ax2.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
ax2.set_title('Azimuthal width', fontsize=6)
ax2.set_xlabel(r'q ($\AA^{-1}$)')
ax2.set_ylabel(r'$\sigma_{\rm azimuth}$ (pixels)')
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.plot(x, y_line, label='fit', color='m', linestyle=':')
ax2.plot(x, y_line_plus_sigmoid, label='fit', color='m')
                        
matplotlib.pyplot.savefig('%s/fits_threshold_%f_ph.png'%(outputFolder, photonThreshold), dpi=4*96)
matplotlib.pyplot.close()