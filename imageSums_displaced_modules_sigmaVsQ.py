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
photonThreshold = 1.0
filter_distanceThreshold = 2.0
nCountsPerPhoton = 26
ellipse_multiplicative_factor = 2.5 
precision_factor = 1
d_threshold = 3

# TO BE FITTED
qs = []
sigXs= []
sigYs= []

Is = []

# FOLDERS
outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber
if not os.path.exists('%s'%outputFolder):
    os.mkdir('%s'%outputFolder)
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
fOpen = open('./Output_r%s/Output_imageSums_moduleDisplacements/imageSums_displaced_modules.txt'%runNumber, 'r')
for hk_line in fOpen:
    splittedLine = hk_line.split()
    try:
        h = int(splittedLine[0])
    except:
        continue
    
    h              = int(splittedLine[0])
    k              = int(splittedLine[1])
    gaussAmplitude = float(splittedLine[2])
    x0             = float(splittedLine[3])
    y0             = float(splittedLine[4])    
    sigX           = float(splittedLine[5])
    sigY           = float(splittedLine[6])
    Igauss         = float(splittedLine[7])
    Isum           = float(splittedLine[8])
    
    reciprocalVector = [h, k]*reciprocalCellRows
    q_x = reciprocalVector[0,0]         # A^(-1)
    q_y = reciprocalVector[0,1]         # A^(-1)
    q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
    
    if gaussAmplitude < 0:
        print '(%4d, %4d) DISCARDED: negative Gauss amplitude'%(h, k)
        continue    
    if abs(Igauss) < photonThreshold*nCountsPerPhoton:
        print '(%4d, %4d) DISCARDED: I < %f photons'%(h, k, photonThreshold)
        continue
    if Isum < photonThreshold:
        print '(%4d, %4d) DISCARDED: I < %f photons'%(h, k, photonThreshold)
        continue
    if sigX > 5 or sigY > 5:
        print '(%4d, %4d) DISCARDED: Gauss width is too large.'%(h, k)
        continue
    if q < 0.15:
        continue
        
    
    if q < 0.9 and sigY > 1.6:
        print '***************', h, k, q, sigY
        
    if sigX < 0:
        print 'NEGATIVE sigma_x', h, k, q, sigX
        continue
    if sigY < 0:
        print 'NEGATIVE sigma_y', h, k, q, sigY
        continue
    

    qs.append(q)
    sigXs.append(sigX)
    sigYs.append(sigY)
    Is.append(Isum)
    
### PLOT SIGMA VS Q ###
matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.savefig('%s/sigXvsQ_%.1fph.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.savefig('%s/sigYvsQ_%.1fph.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()



# FIT sigX VS q WITH A 4th ORDER POLY
popt_sigX, pcov = scipy.optimize.curve_fit(imageSums_utilities.poly_4, qs, sigXs)
x = numpy.linspace(0.1, 1.6, 100)
y_poly_4 = imageSums_utilities.poly_4(x, *popt_sigX)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.title('%s'%popt_sigX)
matplotlib.pyplot.plot(x, y_poly_4, label='fit')
matplotlib.pyplot.savefig('%s/sigXvsQ_poly_4_%.1fph.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()



# FIT sigY VS q WITH A 4th ORDER POLY
popt_sigY, pcov = scipy.optimize.curve_fit(imageSums_utilities.poly_4, qs, sigYs)
x = numpy.linspace(0.1, 1.6, 100)
y_poly_4 = imageSums_utilities.poly_4(x, *popt_sigY)

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.title('%s'%popt_sigY)
matplotlib.pyplot.plot(x, y_poly_4, label='fit')
matplotlib.pyplot.savefig('%s/sigYvsQ_poly_4_%.1fph.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()

offset_sigX = popt_sigX[0]
offset_sigY = popt_sigY[0]

print offset_sigX, offset_sigY

avg_offset = (offset_sigX + offset_sigY)/2
print 'AVERAGE OFFSET: ', avg_offset

# FIT sigX VS q WITH A 4th ORDER POLY, FIXED OFFSET
popt_sigX, pcov = scipy.optimize.curve_fit(lambda (x), a, b: imageSums_utilities.poly_4((x), avg_offset, a, b), (qs), sigXs) 
x = numpy.linspace(0.1, 1.6, 100)
y_poly_4_sigX = imageSums_utilities.poly_4(x, avg_offset, popt_sigX[0], popt_sigX[1])

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.title('%s'%popt_sigX)
matplotlib.pyplot.plot(x, y_poly_4_sigX, label='fit')
matplotlib.pyplot.savefig('%s/sigXvsQ_poly_4_%.1fph_fixed_offset.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()

# FIT sigY VS q WITH A 4th ORDER POLY, FIXED OFFSET
popt_sigY, pcov = scipy.optimize.curve_fit(lambda (x), a, b: imageSums_utilities.poly_4((x), avg_offset, a, b), (qs), sigYs) 
x = numpy.linspace(0.1, 1.6, 100)
y_poly_4_sigY = imageSums_utilities.poly_4(x, avg_offset, popt_sigY[0], popt_sigY[1])

matplotlib.pyplot.figure(figsize = (20, 5))
matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
axes = matplotlib.pyplot.gca()
axes.set_xlim([0, 1.8])
axes.set_ylim([0.5, 4.0])
matplotlib.pyplot.title('%s'%popt_sigY)
matplotlib.pyplot.plot(x, y_poly_4_sigY, label='fit')
matplotlib.pyplot.savefig('%s/sigYvsQ_poly_4_%.1fph_fixed_offset.png'%(outputFolder, photonThreshold))
matplotlib.pyplot.close()



### SUMMARY PLOT ###
myFigure, (ax1, ax2) = matplotlib.pyplot.subplots(2, 1)    
ax1.scatter(qs, sigXs  , color='b', label='data', alpha=0.5, s=3)
#ax1.set_title('Radial width (%s)'%popt_sigX, fontsize=8)
ax1.set_ylabel(r'$\sigma_{\rm radial}$ (pixels)')
ax1.set_ylim([0.5, 4.0])
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.plot(x, y_poly_4_sigX, label='fit', color='m')

ax2.scatter(qs, sigYs, color='b', label='data', alpha=0.5, s=3)
#ax2.set_title('Azimuthal width (%s)'%popt_sigY, fontsize=8)
ax2.set_xlabel(r'q ($\AA^{-1}$)')
ax2.set_ylabel(r'$\sigma_{\rm azimuth}$ (pixels)')
ax2.set_ylim([0.5, 4.0])
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.plot(x, y_poly_4_sigY, label='fit', color='m')

matplotlib.pyplot.tight_layout()                        
matplotlib.pyplot.savefig('%s/fits_poly_4_fixed_offset_threshold_%.1f_ph_paper.png'%(outputFolder, photonThreshold), dpi=4*96)
matplotlib.pyplot.savefig('%s/fits_poly_4_fixed_offset_threshold_%.1f_ph_paper.pdf'%(outputFolder, photonThreshold), dpi=4*96)
matplotlib.pyplot.close()

### REPEAT GAUSS FIT OF IMAGE SUMS USING FIXED SIGMAS ###
fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements/orbits.pkl'%runNumber
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# CHECK HOW MANY ORBITS ARE THERE (220)
nOrbits = 0
for key, orbit in orbitsDictionary.items():
    nOrbits = nOrbits + 1    
print nOrbits, ' ORBITS'

fLog = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse_x0_y0.txt'%outputFolder, 'w')

for key, orbit in orbitsDictionary.items():
    
    label = orbit.label
    h_label = label[0]
    k_label = label[1]
    bgSubtracted_total_sum = orbit.bgSubtracted_total_sum
    q = 2 * numpy.pi / orbit.resolution
    sigma_x = imageSums_utilities.poly_4(q, avg_offset, popt_sigX[0], popt_sigX[1])
    sigma_y = imageSums_utilities.poly_4(q, avg_offset, popt_sigY[0], popt_sigY[1])
    print q, sigma_x, sigma_y
    
    ### SUM INTEGRATION ###
    integratedIntensity_circle = imageSums_utilities.integrate(bgSubtracted_total_sum, expansion_factor=precision_factor)     
    integratedIntensity_circle = integratedIntensity_circle/((10**precision_factor)**2)
    
    integratedIntensity_ellipse = imageSums_utilities.integrate_ellipse(bgSubtracted_total_sum, sigma_x, sigma_y, 
                                                                        ellipse_multiplicative_factor, expansion_factor=precision_factor)
    integratedIntensity_ellipse = integratedIntensity_ellipse/((10**precision_factor)**2)
    
    ### GAUSS FIT AND INTEGRAL ###    
    refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                sigma_x*10**precision_factor, 
                                                                                                                                sigma_y*10**precision_factor)
    if not numpy.isnan(gauss_integral):
        gauss_integral = gauss_integral/((10**precision_factor)**2)
        refined_x0 = refined_x0/(10**precision_factor)
        refined_y0 = refined_y0/(10**precision_factor)
    
    print h_label, k_label, integratedIntensity_circle, integratedIntensity_ellipse, gauss_integral/nCountsPerPhoton 

    ### PLOT GAUSS FIT ###
    if not numpy.isnan(gauss_integral):
        myFigureObject = matplotlib.pyplot.figure()
        myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), 
                                                     origin='lower', interpolation='nearest')
        matplotlib.pyplot.gca().contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                        numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                        data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), 4, colors='w')
        matplotlib.pyplot.title('Orbit: %d %d Gauss integral (fixed sigmas): %.1f counts'
                                 %(h_label, k_label, gauss_integral))
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        
        if orbit.resolution > 7.0:
            matplotlib.pyplot.savefig('%s/low_res/gauss_module_sum_fixed_sigmas_%d_%d.png'
                                       %(outputFolder, h_label, k_label), dpi = 2*96)
        if 6.0 < orbit.resolution <= 7.0:
            matplotlib.pyplot.savefig('%s/high_res_7_6/gauss_module_sum_fixed_sigmas_%d_%d.png'
                                       %(outputFolder, h_label, k_label), dpi = 2*96)
        if 5.0 < orbit.resolution <= 6.0:
            matplotlib.pyplot.savefig('%s/high_res_6_5/gauss_module_sum_fixed_sigmas_%d_%d.png'
                                       %(outputFolder, h_label, k_label), dpi = 2*96)
        if 4.0 < orbit.resolution <= 5.0:
            matplotlib.pyplot.savefig('%s/high_res_5_4/gauss_module_sum_fixed_sigmas_%d_%d.png'
                                       %(outputFolder, h_label, k_label), dpi = 2*96)
        matplotlib.pyplot.close()  

    fLog.write('%4d%4d%8.2f%8.2f%8.2f%8.2f%8.2f\n'%(h_label, k_label, 
                                                    integratedIntensity_circle, gauss_integral/nCountsPerPhoton, integratedIntensity_ellipse, 
                                                    refined_x0, refined_y0))    
fLog.close()

# COMPARE: I_ELLIPSE, I_GAUSS_FIXED_SIGMAS in RESOLUTION RANGE 55 - 4 A (AFTER FILTERING OUT BAD GAUSS FITS)
x0_center = bgSubtracted_total_sum.shape[1]/(2*10**precision_factor)
y0_center = bgSubtracted_total_sum.shape[0]/(2*10**precision_factor)
print 'CENTER', x0_center, y0_center


for I_threshold in [1.0, 0.5, 0.0]:
    Is_sum_circle_python = []
    Is_gauss_fixed_sigmas_python = []
    Is_sum_ellipse_python = []
    
    intensityFile_sum = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse_x0_y0.txt'%outputFolder, 'r')                                  # 55 - 4 A
    for I_line in intensityFile_sum:   # 55 - 4 A
        splittedLine = I_line.split()
        h = int(splittedLine[0])
        k = int(splittedLine[1])
        I_sum_circle_python         = float(splittedLine[2])
        I_gauss_fixed_sigmas_python = float(splittedLine[3])
        I_sum_ellipse_python        = float(splittedLine[4])
        x0                          = float(splittedLine[5])
        y0                          = float(splittedLine[6])
        d = numpy.sqrt((x0-x0_center)**2 + (y0-y0_center)**2)
        
        
        
        if I_sum_circle_python > I_threshold and I_gauss_fixed_sigmas_python > I_threshold and I_sum_ellipse_python > I_threshold and d < d_threshold:
            Is_sum_circle_python.append(I_sum_circle_python)
            Is_gauss_fixed_sigmas_python.append(I_gauss_fixed_sigmas_python)
            Is_sum_ellipse_python.append(I_sum_ellipse_python)
        else:
            print '%d, %d   I < %.1fph'%(h, k, I_threshold)
    intensityFile_sum.close()
                                
    print '\n***** Resolution range 55 A - 4 A  I_thresh = %.2f*****\n'%I_threshold

    CC = imageSums_utilities.Correlate(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
    
    # LINE FIT
    popt_sigY_line, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
    linFit_offset = popt_sigY_line[0]
    linFit_slope  = popt_sigY_line[1]
    print 'Line fit: I_gauss = %.2f + %.2f * I_ellipse'%(linFit_offset, linFit_slope)
    x = numpy.linspace(0.1, 1000, 1000)
    y = x

    print 'I sum (ellipse) python - I gauss (fixed sigmas) python  (n pairs = %d): CC = %.8f'%(len(Is_sum_ellipse_python), CC)
    matplotlib.pyplot.figure()
    matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python), CC))
    matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
    matplotlib.pyplot.plot(x, y, '--', color='r')
    matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - integral on ellipse')
    matplotlib.pyplot.gca().set_ylabel('photon counts - image sum - gauss integral - fixed sigmas', rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_thr_%.1fph.png'%(outputFolder, I_threshold), dpi=3*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()
    #matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse_python), CC))
    
    matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
    matplotlib.pyplot.plot(x, y, '--', color='r')
    matplotlib.pyplot.gca().set_xscale('log')
    matplotlib.pyplot.gca().set_yscale('log')
    matplotlib.pyplot.gca().set_xlim([0.1,1000])
    matplotlib.pyplot.gca().set_ylim([0.1,1000])
    matplotlib.pyplot.axes().set_aspect('equal')
    matplotlib.pyplot.gca().set_xlabel(r'I$_{\rm ellipse}$ (photon counts)')
    matplotlib.pyplot.gca().set_ylabel(r'I$_{\rm gaussian}$ (photon counts)', rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_thr_%.1fph_logscale_CC_%.4f_nPairs_%d.png'%(outputFolder, I_threshold, CC, len(Is_sum_ellipse_python)), dpi=3*96)
    matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_thr_%.1fph_logscale_CC_%.4f_nPairs_%d.pdf'%(outputFolder, I_threshold, CC, len(Is_sum_ellipse_python)), dpi=3*96)
    matplotlib.pyplot.close()

# SAVE SIGMA FIT PARAMETERS
sigmaX_curve = open('%s/sigmaXCurveParameters.pkl'%outputFolder, 'wb')
sigX_parameters = [avg_offset, popt_sigX[0], popt_sigX[1]]
pickle.dump(sigX_parameters, sigmaX_curve)
sigmaX_curve.close()

sigmaY_curve = open('%s/sigmaYCurveParameters.pkl'%outputFolder, 'wb')
sigY_parameters = [avg_offset, popt_sigY[0], popt_sigY[1]]
pickle.dump(sigY_parameters, sigmaY_curve)
sigmaY_curve.close()