# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import pickle
import numpy
import os

import imageSums_utilities

runNumber = '0127'

inputFolder = './Output_r%s/Output_imageSums_moduleDisplacements'%runNumber
outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements/Figures'%runNumber
if not os.path.exists('%s'%outputFolder):
    os.mkdir('%s'%outputFolder)

fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaXCurveParameters.pkl'%runNumber, 'rb')
sigmaXCurveParameters = pickle.load(fRead)                                               
fRead.close()

fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaYCurveParameters.pkl'%runNumber, 'rb')
sigmaYCurveParameters = pickle.load(fRead)                                               
fRead.close()   
    
fileToOpen = '%s/orbits.pkl'%inputFolder
fRead = open(fileToOpen, 'rb')
orbitsDictionary = pickle.load(fRead)                                               
fRead.close()

# IMAGE SUMS SUBPLOTS 
print 'Plot superorbit (2, 11).\n'
orbitIndices = [[-2, -11], [11, 2], [-11, -2], [2, 11]]       
alpha_label = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)']    
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
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[5 : 25, 5 : 25]
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
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[5 : 25, 5 : 25]
                q = 2 * numpy.pi / orbit.resolution
                
                sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
                sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2]) 
                
                refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                    sigma_x, 
                                                                                                                                    sigma_y)
                
    
                ax.set_title('%s'%(alpha_label[idx]), fontsize=8, loc='left')
                ax.set_title('{(%d, %d)}'%(h, k), fontsize=8)
                ax.imshow(data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), origin='lower', interpolation='nearest', vmin=0, vmax=3)
                ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                            numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                            data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1]), 4, colors='w')
                ax.axis('off')
            
                print 'Sigma_x: %.2f, Sigma_y: %.2f'%(sigma_x, sigma_y)

matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11.png'%(outputFolder), dpi = 4*96 )       
matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11.pdf'%(outputFolder), dpi = 4*96 )      
matplotlib.pyplot.close()  

# IMAGE SUMS SUBPLOTS 
print 'Plot superorbit (2, 11).\n'
orbitIndices = [[-2, -11], [11, 2], [-11, -2], [2, 11]]       
alpha_label = ['(f)', '(g)', '(h)', '(i)']    
myFigure, a = matplotlib.pyplot.subplots(1, 4)    
a = a.ravel()
for idx, ax in enumerate(a):   
   
    orbit_label = orbitIndices[idx]
    h = orbit_label[0]
    k = orbit_label[1]
    print orbit_label
    
    for key, orbit in orbitsDictionary.items():
    
        label = orbit.label
        h_label = label[0]
        k_label = label[1]
        
        if h_label == h and k_label == k:
            bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[5 : 25, 5 : 25]
            print bgSubtracted_total_sum.shape
            
            ax.set_title('%s'%(alpha_label[idx]), fontsize=11, loc='left')
            ax.set_title('{(%d, %d)}'%(h, k), fontsize=11)
            ax.imshow(bgSubtracted_total_sum, origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.axis('off')
    
matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i.png'%(outputFolder), dpi = 4*96 )       
matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i.pdf'%(outputFolder), dpi = 4*96 )      
matplotlib.pyplot.close()  



##
print 'Plot Gauss fit and integration area.\n'

orbit_label = [2, 11]
print 'Plotting orbit: ', orbit_label
h = orbit_label[0]
k = orbit_label[1]

for key, orbit in orbitsDictionary.items():
        
    label = orbit.label
    h_label = label[0]
    k_label = label[1]
    
    if h_label == h and k_label == k:
        
        bgSubtracted_total_sum = orbit.bgSubtracted_total_sum
        print 'SHAPE: ', bgSubtracted_total_sum.shape    
        
        q = 2 * numpy.pi / orbit.resolution        
        sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2]) 
        
        refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                    sigma_x, 
                                                                                                                                    sigma_y)
        
        data_toFit = data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
        gaussFit   = data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
        
        # SUMMED IMAGE WITH GAUSS FIT AND ELLIPTICAL INTEGRATION AREA
        matplotlib.pyplot.figure()
        ax = matplotlib.pyplot.gca()        
        ax.imshow(data_toFit, origin='lower', interpolation='nearest', vmin=0, vmax=3)
        ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                   numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                   gaussFit, 
                   4, linestyles='dashed', colors='w', linewidth=3)
        ax.axis('off')
        
        matplotlib.pyplot.plot((refined_x0, refined_x0), (0-0.5, bgSubtracted_total_sum.shape[0]-0.5), color='red', linestyle='dashed', linewidth=2)
        matplotlib.pyplot.plot((0-0.5, bgSubtracted_total_sum.shape[1]-0.5), (refined_y0, refined_y0), color='red', linestyle='dashed', linewidth=2)
                
        mean    = [bgSubtracted_total_sum.shape[1]/2, bgSubtracted_total_sum.shape[0]/2]        
        width   = 2* (2.5*sigma_x)
        height  = 2* (2.5*sigma_y)
        angle   = 0
        ellipse = matplotlib.patches.Ellipse(xy=mean, width=width, height=height, angle=angle, ec='m', fc='none', linewidth=3)
        ax.add_patch(ellipse)
    
        matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d_Final.png'%(outputFolder, h, k), dpi = 4*96 )       
        matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d_Final.pdf'%(outputFolder, h, k), dpi = 4*96 )      
        matplotlib.pyplot.close()  
        
        print 'Gaussian center: ', refined_x0, refined_y0
        print 'Ellipse center: ', mean
        
        # HORIZONTAL EXTRACT (DATA AND FIT)        
        horizontal_extract = data_toFit[int(bgSubtracted_total_sum.shape[0]/2), :]
        ns = numpy.linspace(0, len(horizontal_extract)-1, len(horizontal_extract))
        matplotlib.pyplot.figure()
        matplotlib.pyplot.scatter(ns, horizontal_extract, edgecolor='none')
        horizontal_extract_fit = gaussFit[int(bgSubtracted_total_sum.shape[0]/2), :]
        matplotlib.pyplot.plot(horizontal_extract_fit, 'c', linewidth=2)
        matplotlib.pyplot.gca().set_xlim([-0.5, 29.5])
        matplotlib.pyplot.gca().set_xticklabels([])
        matplotlib.pyplot.gca().set_ylabel(r"$I$ (photons)", fontsize = 22, rotation = 'vertical')
        matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=10)   
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig('%s/Gauss_Fit_horizontalExtract_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
        matplotlib.pyplot.savefig('%s/Gauss_Fit_horizontalExtract_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
        matplotlib.pyplot.close()  
        
        # VERTICAL EXTRACT (DATA AND FIT)
        vertical_extract = data_toFit[:, int(bgSubtracted_total_sum.shape[1]/2)]
        ns = numpy.linspace(0, len(vertical_extract)-1, len(vertical_extract))
        matplotlib.pyplot.figure()
        matplotlib.pyplot.scatter(ns, vertical_extract, edgecolor='none')
        vertical_extract_fit = gaussFit[:, int(bgSubtracted_total_sum.shape[1]/2)]
        matplotlib.pyplot.plot(vertical_extract_fit, 'c', linewidth=2)
        matplotlib.pyplot.gca().set_xlim([-0.5, 29.5])
        matplotlib.pyplot.gca().set_xticklabels([])
        matplotlib.pyplot.gca().set_ylabel(r"$I$ (photons)", fontsize = 22, rotation = 'vertical')
        matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=10)    
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig('%s/Gauss_Fit_verticalExtract_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
        matplotlib.pyplot.savefig('%s/Gauss_Fit_verticalExtract_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
        matplotlib.pyplot.close()

print 'Plotting all orbits.'

for key, orbit in orbitsDictionary.items():
        
    label = orbit.label
    h = label[0]
    k = label[1]
    print h, k 
            
    bgSubtracted_total_sum = orbit.bgSubtracted_total_sum
        
    q = 2 * numpy.pi / orbit.resolution        
    sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
    sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2]) 
    
    refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                sigma_x, 
                                                                                                                                sigma_y)
    
    data_toFit = data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
    gaussFit   = data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
    
    # SUMMED IMAGE WITH GAUSS FIT AND ELLIPTICAL INTEGRATION AREA
    matplotlib.pyplot.figure()
    ax = matplotlib.pyplot.gca()        
    ax.imshow(data_toFit, origin='lower', interpolation='nearest', vmin=0, vmax=3)
    ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
               numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
               gaussFit, 
               4, linestyles='dashed', colors='w', linewidth=3)
    ax.axis('off')
    matplotlib.pyplot.plot((refined_x0, refined_x0), (0-0.5, bgSubtracted_total_sum.shape[0]-0.5), color='red', linestyle='dashed', linewidth=2)
    matplotlib.pyplot.plot((0-0.5, bgSubtracted_total_sum.shape[1]-0.5), (refined_y0, refined_y0), color='red', linestyle='dashed', linewidth=2)
         
    mean    = [bgSubtracted_total_sum.shape[1]/2, bgSubtracted_total_sum.shape[0]/2]
    width   = 2* (2.5*sigma_x)
    height  = 2* (2.5*sigma_y)
    angle   = 0
    ellipse = matplotlib.patches.Ellipse(xy=mean, width=width, height=height, angle=angle, ec='m', fc='none', linewidth=3)
    ax.add_patch(ellipse)

    matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
    matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
    matplotlib.pyplot.close()  
    print 'Gaussian center: ', refined_x0, refined_y0
    print 'Ellipse center: ', mean