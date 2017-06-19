# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import pickle
import numpy

import imageSums_utilities

runNumber = '0127'
precision_factor = 1

outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements'%runNumber

fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaXCurveParameters.pkl'%runNumber, 'rb')
sigmaXCurveParameters = pickle.load(fRead)                                               
fRead.close()

fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaYCurveParameters.pkl'%runNumber, 'rb')
sigmaYCurveParameters = pickle.load(fRead)                                               
fRead.close()   


orbitIndices = [[-2, -11], [11, 2], [-11, -2], [2, 11]]       
alpha_label = ['(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)'] 

    
fileToOpen = '%s/orbits.pkl'%outputFolder
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# IMAGE SUMS SUBPLOTS    
myFigure, a = matplotlib.pyplot.subplots(2, 4)    
a = a.ravel()
for idx, ax in enumerate(a):   
    
    index = idx%4
    
    orbit_label = orbitIndices[index]
    h = orbit_label[0]
    k = orbit_label[1]
    print orbit_label
    if idx < 4:
        for key, orbit in orbitsDictionary.items():
        
            label = orbit.label
            h_label = label[0]
            k_label = label[1]
            
            if h_label == h and k_label == k:
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
                print bgSubtracted_total_sum.shape
                
                ax.set_title('%s'%(alpha_label[idx]), fontsize=8, loc='left')
                ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
                ax.imshow(bgSubtracted_total_sum, origin='lower', interpolation='nearest', vmin=0, vmax=3)
                ax.axis('off')
    else:
        for key, orbit in orbitsDictionary.items():
        
            label = orbit.label
            h_label = label[0]
            k_label = label[1]
            
            if h_label == h and k_label == k:
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
                q = 2 * numpy.pi / orbit.resolution
                
                sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
                sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2]) 
                
                refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                    sigma_x*10**precision_factor, 
                                                                                                                                    sigma_y*10**precision_factor)
                
    
                ax.set_title('%s'%(alpha_label[idx]), fontsize=8, loc='left')
                ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
                ax.imshow(data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), origin='lower', interpolation='nearest', vmin=0, vmax=3)
                ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                            numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                            data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), 4, colors='w')
                ax.axis('off')
            
                print 'Sigma_x: %.2f, Sigma_y: %.2f'%(sigma_x, sigma_y)

matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_f_to_m.png'%(outputFolder), dpi = 4*96 )       
matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_f_to_m.pdf'%(outputFolder), dpi = 4*96 )      
matplotlib.pyplot.close()  