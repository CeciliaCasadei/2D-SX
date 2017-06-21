# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import scipy.optimize
import pickle

import imageSums_utilities



# PARAMETERS
runNumber = '0127'
precision_factor = 1
d_threshold = 3

# FOLDERS
outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber

## LOAD DATA
fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements/orbits.pkl'%runNumber
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# CHECK HOW MANY ORBITS ARE THERE (220)
nOrbits = 0
for key, orbit in orbitsDictionary.items():
    bgSubtracted_total_sum = orbit.bgSubtracted_total_sum
    nOrbits = nOrbits + 1    
print nOrbits, ' ORBITS'

# GET CENTER COORDINATES
x0_center = bgSubtracted_total_sum.shape[1]/(2*10**precision_factor)
y0_center = bgSubtracted_total_sum.shape[0]/(2*10**precision_factor)
print 'CENTER', x0_center, y0_center

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