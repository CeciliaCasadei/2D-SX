# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import os
import numpy
import scipy.optimize
import pickle

import imageSums_utilities

# PARAMETERS
photonThreshold = 2.0
nCountsPerPhoton = 26

# TO BE FITTED
qs = []
sigXs= []
sigYs= []

Is = []

# FOLDERS
outputFolder = './Output_imageSums_sigmaFits/paper_figure'
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
        print '(%4d, %4d) DISCARDED: I < %.2f photons'%(h, k, photonThreshold)
        continue
    if sigX > 7 or sigY > 7:
        print '(%4d, %4d) DISCARDED: Gauss width is too large.'%(h, k)
        continue

    qs.append(q)
    sigXs.append(sigX)
    sigYs.append(sigY)
    Is.append(Isum)

# FIT sigX VS q WITH A QUADRATIC
popt_sigX, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic, qs, sigXs)
x = numpy.linspace(0.1, 1.6, 100)
y_quadratic = imageSums_utilities.quadratic(x, *popt_sigX)

# FIT sigY VS q WITH A LINE AT LOW q
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.9]
low_q_sigYs = [sigYs[i] for i in range(0, len(qs)) if qs[i] < 0.9]

popt_sigY_line, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, low_qs, low_q_sigYs)
y_line = imageSums_utilities.line(x, *popt_sigY_line)

offset = popt_sigY_line[0]
slope = popt_sigY_line[1]

# FIT sigY VS q WITH A LINE + SIGMA
popt_sigY, pcov = scipy.optimize.curve_fit(lambda x, x0, k, scale: imageSums_utilities.line_plus_sigmoid(x, offset, slope, x0, k, scale), qs, sigYs)
y_line_plus_sigmoid = imageSums_utilities.line_plus_sigmoid(x, offset, slope, popt_sigY[0], popt_sigY[1], popt_sigY[2])

### 2x1 PLOT ###
myFigure, (ax1, ax2) = matplotlib.pyplot.subplots(2, 1)    
ax1.scatter(qs, sigXs  , color='b', label='data', alpha=0.5)
ax1.set_title('Radial width', fontsize=8)
ax1.set_ylabel(r'$\sigma_{\rm radial}$ (pixels)')
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.plot(x, y_quadratic, label='fit', color='m')

ax2.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
ax2.set_title('Azimuthal width', fontsize=8)
ax2.set_xlabel(r'q ($\AA^{-1}$)')
ax2.set_ylabel(r'$\sigma_{\rm azimuth}$ (pixels)')
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.plot(x, y_line, label='fit', color='m', linestyle=':')
ax2.plot(x, y_line_plus_sigmoid, label='fit', color='m')
                        
matplotlib.pyplot.savefig('%s/fits_threshold_%f_ph.png'%(outputFolder, photonThreshold), dpi=4*96)
matplotlib.pyplot.close()

### PAPER FIGURE ###

# DEFINE ORBIT INDICES 
orbitIndices = [[2, 11], [11, 2], [-2, -11], [-11, -2]]       
alpha_label = ['(i)', '(j)', '(k)', '(l)'] 

# GAUSS FIT USING FIXED SIGMAS
fileToOpen = './Output_imageSums/orbits.pkl'
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

myFigure, a = matplotlib.pyplot.subplots(1, 4)    
a = a.ravel()
for idx, ax in enumerate(a):   
    print idx
    orbit_label = orbitIndices[idx]
    h = orbit_label[0]
    k = orbit_label[1]
    for key, orbit in orbitsDictionary.items():
    
        label = orbit.label
        h_label = label[0]
        k_label = label[1]
        
        if h_label == h and k_label == k:
            bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtractedOrbitSumMatrix_rotated[5:25, 5:25]
            q = 2 * numpy.pi / orbit.resolution
            sigma_x = imageSums_utilities.quadratic(q, *popt_sigX)
            sigma_y = imageSums_utilities.line_plus_sigmoid(q, offset, slope, popt_sigY[0], popt_sigY[1], popt_sigY[2])  
                     
            n_x = bgSubtractedOrbitSumMatrix_rotated.shape[1]
            n_y = bgSubtractedOrbitSumMatrix_rotated.shape[0]
            x = numpy.linspace(0, n_x-1, n_x)
            y = numpy.linspace(0, n_y-1, n_y)
            x, y = numpy.meshgrid(x, y)                      # column_index, row_index
            
            data = bgSubtractedOrbitSumMatrix_rotated.ravel()           
            data = data.T
            data = numpy.asarray(data)
            data = data.flatten()
            
            initial_x = float(n_x)/2
            initial_y = float(n_y)/2
            initial_guess = (numpy.amax(bgSubtractedOrbitSumMatrix_rotated), initial_x, initial_y)
                    
            popt, pcov = scipy.optimize.curve_fit(lambda (x, y), amplitude, xo, yo: imageSums_utilities.twoD_Gaussian_simple((x, y), amplitude, xo, yo, sigma_x, sigma_y), (x, y), data, p0=initial_guess)
           
            refined_amplitude = popt[0]
            refined_x0 = popt[1]
            refined_y0 = popt[2]
    
            data_fitted = imageSums_utilities.twoD_Gaussian_simple((x, y), refined_amplitude, refined_x0, refined_y0, sigma_x, sigma_y) 
             
            ax.set_title('%s'%(alpha_label[idx]), fontsize=8, loc='left')
            ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
            ax.imshow(data.reshape(n_y, n_x), origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.contour(x, y, data_fitted.reshape(n_y, n_x), 4, colors='w')
            ax.axis('off')
            
            print 'Sigma_x: %.2f, Sigma_y: %.2f'%(sigma_x, sigma_y)
              
matplotlib.pyplot.savefig('%s/paper_x4_bgsub_rotatetd_gauss_fixed_sigma_orbit_2_11.png'%(outputFolder), dpi = 2*96 )       
matplotlib.pyplot.close()  

