# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import scipy.optimize
import warnings

import imageSums_utilities


def plot(myArguments):
    str1 = "--selectedRun <selectedRun> --d_threshold <d_threshold> --halfWidth <halfWidth>"
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "d_threshold=",
                                                                 "halfWidth="])
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules_integrationPlots.py %s'%str1
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules_integrationPlots.py %s'%str1
            sys.exit()
        elif option == "--selectedRun":
            runNumber = value.zfill(4)
        elif option == "--d_threshold":
            d_threshold = float(value)
        elif option == "--halfWidth":
            halfWidth = int(value)
    
    warnings.filterwarnings("ignore")
    
    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber
        
    # CENTER COORDINATES
    x0_center = halfWidth-10
    y0_center = halfWidth-10
    
    # COMPARE: I_ELLIPSE, I_GAUSS_FIXED_SIGMAS in RESOLUTION RANGE 55 - 4 A (AFTER FILTERING OUT BAD GAUSS FITS)
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
            d                           = numpy.sqrt((x0-x0_center)**2 + (y0-y0_center)**2)
            
            if I_sum_circle_python > I_threshold and I_gauss_fixed_sigmas_python > I_threshold and I_sum_ellipse_python > I_threshold and d < d_threshold:
                Is_sum_circle_python.append(I_sum_circle_python)
                Is_gauss_fixed_sigmas_python.append(I_gauss_fixed_sigmas_python)
                Is_sum_ellipse_python.append(I_sum_ellipse_python)
            else:
                if d >= d_threshold:
                    print '%6d, %6d   d > %.1f pxls'%(h, k, d_threshold)
                if I_sum_circle_python <= I_threshold or I_gauss_fixed_sigmas_python <= I_threshold or I_sum_ellipse_python <= I_threshold:
                    print '%6d, %6d   I < %.1f ph'%(h, k, I_threshold)
                    
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
        matplotlib.pyplot.scatter(Is_sum_ellipse_python, Is_gauss_fixed_sigmas_python)
        matplotlib.pyplot.plot(x, y, '--', color='r')
        matplotlib.pyplot.gca().set_xscale('log')
        matplotlib.pyplot.gca().set_yscale('log')
        matplotlib.pyplot.gca().set_xlim([I_threshold-0.3,500])
        matplotlib.pyplot.gca().set_ylim([I_threshold-0.3,500])
        matplotlib.pyplot.axes().set_aspect('equal')
        matplotlib.pyplot.gca().set_xlabel(r'I$_{\rm ellipse}$ (photon counts)')
        matplotlib.pyplot.gca().set_ylabel(r'I$_{\rm gaussian}$ (photon counts)', rotation = 'vertical')
        matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_thr_%.1fph_logscale_CC_%.4f_nPairs_%d.png'%(outputFolder, I_threshold, CC, len(Is_sum_ellipse_python)), dpi=3*96)
        matplotlib.pyplot.savefig('%s/Is_gauss_fixed_sigmas_python_VS_Is_sum_ellipse_python_4A_thr_%.1fph_logscale_CC_%.4f_nPairs_%d.pdf'%(outputFolder, I_threshold, CC, len(Is_sum_ellipse_python)), dpi=3*96)
        matplotlib.pyplot.close()
    
    
    
if __name__ == "__main__":
    print "\n**** CALLING imageSums_displaced_modules_integrationPlots ****"
    plot(sys.argv[1:])