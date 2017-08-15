# -*- coding: utf-8 -*-
### REPEAT GAUSS FIT OF IMAGE SUMS USING FIXED SIGMAS ###


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle
import sys
import getopt

import imageSums_utilities

def finalIntegration(myArguments):
    
    str1 = "--selectedRun <selectedRun> --nCountsPerPhoton <nCountsPerPhoton> --ellipse_multiplicative_factor <ellipse_multiplicative_factor> --precisionFactor <precisionFactor>"
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "nCountsPerPhoton=", 
                                                                 "ellipse_multiplicative_factor=",
                                                                 "precisionFactor="])
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules_finalIntegration_pixelConversion.py %s'%str1
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules_finalIntegration_pixelConversion.py %s'%str1
            sys.exit()
        elif option == "--selectedRun":
            runNumber = value.zfill(4)
        elif option == "--nCountsPerPhoton":
            nCountsPerPhoton = int(value)
        elif option == "--ellipse_multiplicative_factor":
            ellipse_multiplicative_factor = float(value)
        elif option == "--precisionFactor":
            precision_factor = int(value)

    
    
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

if __name__ == "__main__":
    print "\n**** CALLING imageSums_displaced_modules_finalIntegration_pixelConversion ****"
    finalIntegration(sys.argv[1:])   