# -*- coding: utf-8 -*-
### REPEAT GAUSS FIT OF IMAGE SUMS USING FIXED SIGMAS ###


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



# FOLDERS
outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber
if not os.path.exists('%s/high_res_7_6'%outputFolder):
        os.mkdir('%s/high_res_7_6'%outputFolder)
if not os.path.exists('%s/high_res_6_5'%outputFolder):
    os.mkdir('%s/high_res_6_5'%outputFolder)
if not os.path.exists('%s/high_res_5_4'%outputFolder):
    os.mkdir('%s/high_res_5_4'%outputFolder)
if not os.path.exists('%s/low_res'%outputFolder):
    os.mkdir('%s/low_res'%outputFolder)



# LOAD DATA
fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements/orbits.pkl'%runNumber
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

fRead = open('%s/sigmaXCurveParameters.pkl'%outputFolder, 'rb')
sigmaXCurveParameters = pickle.load(fRead)                                               
fRead.close()

fRead = open('%s/sigmaYCurveParameters.pkl'%outputFolder, 'rb')
sigmaYCurveParameters = pickle.load(fRead)                                               
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
    sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
    sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
    print q, sigma_x, sigma_y
    
    ### SUM INTEGRATION ON FIXED RADIUS CIRCLE ###
    integratedIntensity_circle = imageSums_utilities.integrate(bgSubtracted_total_sum, expansion_factor=precision_factor)     
    integratedIntensity_circle = integratedIntensity_circle/((10**precision_factor)**2)
    
    ### SUM INTEGRATION ON VARIABLE ELLIPSE ###
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