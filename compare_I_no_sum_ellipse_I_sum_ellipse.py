# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

import imageSums_utilities

# COMPARE: I_NO_SUM_ELLIPSE_PY, I_SUM_ELLIPSE_PY in RESOLUTION RANGE 55 - 4 A
Is_sum_ellipse = []
Is_no_sum_ellipse = []
Is_no_sum_ellipse_error = []

intensityFile_sum_ellipse    = open('./Output_imageSums_sigmaFits/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse.txt', 'r') # 55 - 4 A
intensityFile_no_sum_ellipse = open('./Output_elliptical_integration/mergedIntensities_h_k_qRod_l_avgI_sigmaI.txt', 'r')         # 55 - 4 A

intensityFile_no_sum_ellipse_list = list(intensityFile_no_sum_ellipse)
nList = len(intensityFile_no_sum_ellipse_list)

for I_line in intensityFile_sum_ellipse:   # 55 - 4 A    h k I_sum_circle I_sum_gauss_fixed_sigmas I_sum_ellipse
    splittedLine = I_line.split()
    h = int(splittedLine[0])
    k = int(splittedLine[1])
    I_sum_ellipse = float(splittedLine[4])
    I_sum_gauss = float(splittedLine[3])
    if not numpy.isnan(I_sum_gauss):    
        for i in range(0, nList):        
            listElement = intensityFile_no_sum_ellipse_list[i]        # h k qRod l avgI_no_sum_ellipse sigI
            splitElement = listElement.split()
            if int(splitElement[0]) == h and int(splitElement[1]) == k:
                I_no_sum_ellipse = float(splitElement[4]) 
                I_no_sum_ellipse_error = float(splitElement[5])
                Is_sum_ellipse.append(I_sum_ellipse)
                Is_no_sum_ellipse.append(I_no_sum_ellipse)
                Is_no_sum_ellipse_error.append(I_no_sum_ellipse_error)
                        
print '\n***** Resolution range 55 A - 4 A *****\n' 
CC = imageSums_utilities.Correlate(Is_sum_ellipse, Is_no_sum_ellipse)
print 'I sum ellipse - I no sum, ellipse integration (n pairs = %d): CC = %.4f'%(len(Is_sum_ellipse), CC)

intensityFile_sum_ellipse.close()
intensityFile_no_sum_ellipse.close()
       
matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse), CC), y=1.05)
matplotlib.pyplot.scatter(Is_sum_ellipse, Is_no_sum_ellipse)
matplotlib.pyplot.gca().set_xlabel('photon counts - image sum - ellipse integral')
matplotlib.pyplot.gca().set_ylabel('photon counts - no sum - ellipse integral', rotation = 'vertical')
matplotlib.pyplot.savefig('./Output_elliptical_integration/Is_sum_ellipse_python_VS_Is_no_sum_ellipse_python_4A.png')
matplotlib.pyplot.close()

matplotlib.pyplot.figure()
matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_sum_ellipse), CC), y=1.05)
matplotlib.pyplot.scatter(Is_sum_ellipse, Is_no_sum_ellipse)
matplotlib.pyplot.gca().set_xscale('log')
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_xlim([0.1, 1000])
matplotlib.pyplot.gca().set_ylim([0.1, 1000])
matplotlib.pyplot.gca().set_ylabel(r'I$_{\rm merged, ellipse}$ (photon counts)')
matplotlib.pyplot.gca().set_xlabel(r'I$_{\rm sum, ellipse}$ (photon counts)')
matplotlib.pyplot.axes().set_aspect('equal')
matplotlib.pyplot.savefig('./Output_elliptical_integration/Is_sum_ellipse_python_VS_Is_no_sum_ellipse_python_4A_logscale.png', dpi=96*4)
matplotlib.pyplot.close()
