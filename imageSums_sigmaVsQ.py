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
runNumber = '0127'
photonThreshold = 2.0
filter_distanceThreshold = 2.0
nCountsPerPhoton = 26
ellipse_multiplicative_factor = 2.5

# TO BE FITTED
qs = []
sigXs= []
sigYs= []

Is = []

# FOLDERS
inputFolder = './Output_r%s/Output_imageSums'%runNumber
outputFolder = './Output_r%s/Output_imageSums_sigmaFits'%runNumber
outputFigures = '%s/Figures'%outputFolder
if not os.path.exists('%s'%outputFolder):
    os.mkdir('%s'%outputFolder)
if not os.path.exists('%s'%outputFigures):
    os.mkdir('%s'%outputFigures)
if not os.path.exists('%s/high_res_7_6'%outputFolder):
        os.mkdir('%s/high_res_7_6'%outputFolder)
if not os.path.exists('%s/high_res_6_5'%outputFolder):
    os.mkdir('%s/high_res_6_5'%outputFolder)
if not os.path.exists('%s/high_res_5_4'%outputFolder):
    os.mkdir('%s/high_res_5_4'%outputFolder)
if not os.path.exists('%s/low_res'%outputFolder):
    os.mkdir('%s/low_res'%outputFolder)

# RECIPROCAL CELL
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I                                                  # A^(-1)

# LOOP ON SPOTS WITH GAUSS INTENSITY ABOVE 1 PHOTON AND GOOD GAUSSIAN FIT
fOpen = open('%s/h_k_Isum_Igauss_sigX_sigY.txt'%inputFolder, 'r')
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
    if sigX > 6 or sigY > 6:
        print '(%4d, %4d) DISCARDED: Gauss width is too large.'%(h, k)
        continue

    qs.append(q)
    sigXs.append(sigX)
    sigYs.append(sigY)
    Is.append(Isum)

### FIND SIGMA VS Q BEHAVIOUR ###

# FIT sigX VS q WITH A QUADRATIC
popt_sigX, pcov = scipy.optimize.curve_fit(imageSums_utilities.quadratic, qs, sigXs)
x = numpy.linspace(0.1, 1.6, 100)
y_quadratic = imageSums_utilities.quadratic(x, *popt_sigX)
matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_quadratic, label='fit')
matplotlib.pyplot.savefig('%s/sigXvsQ_quadraticFit.png'%outputFigures)
matplotlib.pyplot.close()

# FIT sigY VS q WITH A LINE AT LOW q
low_qs      = [qs[i]    for i in range(0, len(qs)) if qs[i] < 0.9]
low_q_sigYs = [sigYs[i] for i in range(0, len(qs)) if qs[i] < 0.9]

popt_sigY_line, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, low_qs, low_q_sigYs)
y_line = imageSums_utilities.line(x, *popt_sigY_line)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
matplotlib.pyplot.plot(x, y_line, label='fit')
matplotlib.pyplot.savefig('%s/sigYvsQ_lineFit.png'%outputFigures)
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
matplotlib.pyplot.savefig('%s/sigYvsQ_linePlusSigmidFit.png'%outputFigures)
matplotlib.pyplot.close()

### SUMMARY PLOT ###
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
                        
matplotlib.pyplot.savefig('%s/fits_threshold_%f_ph.png'%(outputFigures, photonThreshold), dpi=4*96)
matplotlib.pyplot.close()

### REPEAT GAUSS FIT OF IMAGE SUMS USING FIXED SIGMAS ###
fileToOpen = '%s/orbits.pkl'%inputFolder
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# CHECK HOW MANY ORBITS ARE THERE (220)
nOrbits = 0
for key, orbit in orbitsDictionary.items():
    nOrbits = nOrbits + 1    
print nOrbits

fLog = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse.txt'%outputFolder, 'w')

for key, orbit in orbitsDictionary.items():
    
    label = orbit.label
    h_label = label[0]
    k_label = label[1]
    bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtractedOrbitSumMatrix_rotated
    q = 2 * numpy.pi / orbit.resolution
    sigma_x = imageSums_utilities.quadratic(q, *popt_sigX)
    sigma_y = imageSums_utilities.line_plus_sigmoid(q, offset, slope, popt_sigY[0], popt_sigY[1], popt_sigY[2])    
    print bgSubtractedOrbitSumMatrix_rotated.shape, ellipse_multiplicative_factor*sigma_x, ellipse_multiplicative_factor*sigma_y
        
    ### SUM INTEGRATION ###
    integratedIntensity_circle  = imageSums_utilities.integrate(bgSubtractedOrbitSumMatrix_rotated)
    integratedIntensity_ellipse = imageSums_utilities.integrate_ellipse(bgSubtractedOrbitSumMatrix_rotated, sigma_x, sigma_y, ellipse_multiplicative_factor)
     
    ### GAUSS FIT AND INTEGRAL ###
    try:
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

        ### ANALYTICAL GAUSSIAN INTEGRAL ###
        gauss_integral = 2 * numpy.pi * refined_amplitude * sigma_x * sigma_y
                    
        ### FILTER OUT BAD FITS ###
        if refined_amplitude < 0:
            print '%4d, %4d Gaussian amplitude < 0'%(h_label, k_label)
            gauss_integral = numpy.nan
        if numpy.sqrt((refined_x0-initial_x)**2 + (refined_y0-initial_y)**2) > filter_distanceThreshold:
            print '%4d, %4d Gaussian center far from prediction'%(h_label, k_label)
            gauss_integral = numpy.nan
        
        ## PLOT GAUSS FIT ###
        myFigureObject = matplotlib.pyplot.figure()
        myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(30, 30), origin='lower', interpolation='nearest')
        data_fitted = imageSums_utilities.twoD_Gaussian_simple((x, y), refined_amplitude, refined_x0, refined_y0, sigma_x, sigma_y) 
        matplotlib.pyplot.gca().contour(x, y, data_fitted.reshape(30, 30), 4, colors='w')
        matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f Gauss integral: %.1f counts'%(h_label, k_label, orbit.resolution, gauss_integral))
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        
        if orbit.resolution > 7.0:
            matplotlib.pyplot.savefig('%s/low_res/gauss_bgsub_rotatetd_orbit_%d_%d.png'%(outputFolder, h_label, k_label), dpi = 2*96)
        if 6.0 < orbit.resolution <= 7.0:
            matplotlib.pyplot.savefig('%s/high_res_7_6/gauss_bgsub_rotated_orbit_%d_%d.png'%(outputFolder, h_label, k_label), dpi = 2*96)
        if 5.0 < orbit.resolution <= 6.0:
            matplotlib.pyplot.savefig('%s/high_res_6_5/gauss_bgsub_rotated_orbit_%d_%d.png'%(outputFolder, h_label, k_label), dpi = 2*96)
        if 4.0 < orbit.resolution <= 5.0:
            matplotlib.pyplot.savefig('%s/high_res_5_4/gauss_bgsub_rotated_orbit_%d_%d.png'%(outputFolder, h_label, k_label), dpi = 2*96)
        matplotlib.pyplot.close()  
        
    except:
        print '%4d, %4d Gaussian fit not possible'%(h_label, k_label)
        gauss_integral = numpy.nan
    
    fLog.write('%4d%4d%8.2f%8.2f%8.2f\n'%(h_label, k_label, integratedIntensity_circle, gauss_integral/nCountsPerPhoton, integratedIntensity_ellipse))    
fLog.close()

# COMPARE: I_SUM_PY, I_GAUSS_FIXED_SIGMAS_PY in RESOLUTION RANGE 55 - 4 A (AFTER FILTERING OUT BAD GAUSS FITS)
Is_sum_circle_python = []
Is_gauss_fixed_sigmas_python = []
Is_sum_ellipse_python = []

intensityFile_sum = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse.txt'%outputFolder, 'r')                                  # 55 - 4 A
for I_line in intensityFile_sum:   # 55 - 4 A
    splittedLine = I_line.split()
    h = int(splittedLine[0])
    k = int(splittedLine[1])
    I_sum_circle_python         = float(splittedLine[2])
    I_gauss_fixed_sigmas_python = float(splittedLine[3])
    I_sum_ellipse_python        = float(splittedLine[4])
    if not numpy.isnan(I_gauss_fixed_sigmas_python):
        Is_sum_circle_python.append(I_sum_circle_python)
        Is_gauss_fixed_sigmas_python.append(I_gauss_fixed_sigmas_python)
        Is_sum_ellipse_python.append(I_sum_ellipse_python)
    else:
        print 'I_gauss_fixed_sigmas_python isnan:', h, k
intensityFile_sum.close()
                            
print '\n***** Resolution range 55 A - 4 A *****\n' 
CC = imageSums_utilities.Correlate(Is_sum_circle_python, Is_gauss_fixed_sigmas_python)
print 'I sum (circle) python  - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_gauss_fixed_sigmas_python), CC)

matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_gauss_fixed_sigmas_python), CC))
matplotlib.pyplot.scatter(Is_sum_circle_python, Is_gauss_fixed_sigmas_python)
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on circle')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_circle_python_4A.png'%outputFigures)
matplotlib.pyplot.close()

CC = imageSums_utilities.Correlate(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
print 'I sum (ellipse) python - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_sum_ellipse_python), CC)
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python), CC))
matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on ellipse')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A.png'%outputFigures)
matplotlib.pyplot.close()

matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python), CC))
matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
matplotlib.pyplot.gca().set_xscale('log')
matplotlib.pyplot.gca().set_yscale('log')
#matplotlib.pyplot.gca().set_xlim([0.01,1000])
#matplotlib.pyplot.gca().set_ylim([0.01,1000])
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on ellipse')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_logscale.png'%outputFigures)
matplotlib.pyplot.close()

# COMPARE: I_SUM_PY, I_GAUSS_FIXED_SIGMAS_PY in RESOLUTION RANGE 55 - 4 A (AFTER FILTERING OUT BAD GAUSS FITS + EXCLUDE PH COUNTS < 1 PH)
Is_sum_circle_python_filtered         = [Is_sum_circle_python[i]         for i in range(0, len(Is_sum_circle_python)) if Is_sum_ellipse_python[i] > 1]
Is_gauss_fixed_sigmas_python_filtered = [Is_gauss_fixed_sigmas_python[i] for i in range(0, len(Is_sum_circle_python)) if Is_sum_ellipse_python[i] > 1]
Is_sum_ellipse_python_filtered        = [Is_sum_ellipse_python[i]        for i in range(0, len(Is_sum_circle_python)) if Is_sum_ellipse_python[i] > 1]

print '\n***** Resolution range 55 A - 4 A  - Photon counts > 1 *****\n' 
CC = imageSums_utilities.Correlate(Is_sum_circle_python_filtered, Is_gauss_fixed_sigmas_python_filtered)
print 'I sum (circle) python  - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_gauss_fixed_sigmas_python_filtered), CC)

matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_gauss_fixed_sigmas_python_filtered), CC))
matplotlib.pyplot.scatter(Is_sum_circle_python_filtered, Is_gauss_fixed_sigmas_python_filtered)
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on circle')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_circle_python_4A_phCountsAboveOne.png'%outputFigures)
matplotlib.pyplot.close()

CC = imageSums_utilities.Correlate(Is_sum_ellipse_python_filtered, Is_gauss_fixed_sigmas_python_filtered)
print 'I sum (ellipse) python - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_sum_ellipse_python_filtered), CC)
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python_filtered), CC))
matplotlib.pyplot.scatter(Is_sum_ellipse_python_filtered, Is_gauss_fixed_sigmas_python_filtered)
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on ellipse')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_phCountsAboveOne.png'%outputFigures)
matplotlib.pyplot.close()

matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python_filtered), CC))
matplotlib.pyplot.scatter(Is_sum_ellipse_python_filtered, Is_gauss_fixed_sigmas_python_filtered)
#matplotlib.pyplot.gca().set_xlim([0.1, 1000])
matplotlib.pyplot.gca().set_xscale('log')
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on ellipse')
matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_logscale_phCountsAboveOne.png'%outputFigures)
matplotlib.pyplot.close()

# SAVE SIGMA FIT PARAMETERS
sigmaX_curve = open('%s/sigmaXCurveParameters.pkl'%outputFolder, 'wb')
pickle.dump(popt_sigX, sigmaX_curve)
sigmaX_curve.close()

sigmaY_curve = open('%s/sigmaYCurveParameters.pkl'%outputFolder, 'wb')
sigY_parameters = [offset, slope, popt_sigY[0], popt_sigY[1], popt_sigY[2]]
pickle.dump(sigY_parameters, sigmaY_curve)
sigmaY_curve.close()