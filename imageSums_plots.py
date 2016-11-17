# -*- coding: utf-8 -*-
import matplotlib.pyplot
import pickle
import numpy
import scipy.optimize

import imageSums_utilities



outputFolder = './Output_imageSums'

# DEFINE ORBIT INDICES 
orbitIndices = [[2, 11], [11, 2], [-2, -11], [-11, -2]]       
alpha_label_1 = ['(a)', '(b)', '(c)', '(d)'] 
alpha_label_2 = ['(e)', '(f)', '(g)', '(h)'] 
    
fileToOpen = '%s/orbits.pkl'%outputFolder
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# IMAGE SUMS SUBPLOTS    
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
            print bgSubtractedOrbitSumMatrix_rotated.shape
            
            ax.set_title('%s'%(alpha_label_1[idx]), fontsize=8, loc='left')
            ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
            ax.imshow(bgSubtractedOrbitSumMatrix_rotated, origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.axis('off')
              
matplotlib.pyplot.savefig('%s/paper_x4_bgsub_rotatetd_orbit_2_11_b.png'%(outputFolder), dpi = 2*96 )       
matplotlib.pyplot.close()  




# IMAGE SUMS WITH GAUSS FIT SUBPLOTS
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
            initial_guess = (numpy.amax(bgSubtractedOrbitSumMatrix_rotated), initial_x, initial_y, 2.0, 2.0)
        
            popt, pcov = scipy.optimize.curve_fit(imageSums_utilities.twoD_Gaussian_simple, (x, y), data, p0=initial_guess)
            data_fitted = imageSums_utilities.twoD_Gaussian_simple((x, y), *popt)       
            
            ax.set_title('%s'%(alpha_label_2[idx]), fontsize=8, loc='left')
            ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
            ax.imshow(data.reshape(n_y, n_x), origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.contour(x, y, data_fitted.reshape(n_y, n_x), 4, colors='w')
            ax.axis('off')
            
            refined_sigma_x = popt[3]
            refined_sigma_y = popt[4]
            print 'Sigma_x: %.2f, Sigma_y: %.2f'%(refined_sigma_x, refined_sigma_y)
                
matplotlib.pyplot.savefig('%s/paper_x4_gauss_bgsub_rotatetd_orbit_2_11_b.png'%(outputFolder), dpi = 2*96 )       
matplotlib.pyplot.close()  