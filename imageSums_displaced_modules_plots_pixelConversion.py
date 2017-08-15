# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import pickle
import numpy
import os
import warnings

import imageSums_utilities

def convert(bgSubtracted_total_sum, precision_factor):
    bgSubtracted_total_sum_converted = numpy.zeros(shape=bgSubtracted_total_sum.shape)
    i_limit = bgSubtracted_total_sum.shape[0]
    j_limit = bgSubtracted_total_sum.shape[1]
    m_limit = i_limit/(10**precision_factor)
    n_limit = j_limit/(10**precision_factor)
    for m in range(0, m_limit):
        for n in range(0, n_limit):
            collection = []
            for i in range(0, i_limit):
                for j in range(0, j_limit):
                    if (int(i)/(10**precision_factor) == m) and (int(j)/(10**precision_factor) == n):
                        element = bgSubtracted_total_sum[i, j]
                        collection.append(element)
            if len(collection) != 100:
                print 'PROBLEM'
            collection = numpy.asarray(collection)
            avg = numpy.average(collection)
            for i in range(0, i_limit):
                for j in range(0, j_limit):
                    if (int(i)/(10**precision_factor) == m) and (int(j)/(10**precision_factor) == n):
                        bgSubtracted_total_sum_converted[i, j] = avg
    return bgSubtracted_total_sum_converted

def plot(myArguments):
    
    str1 = "--selectedRun <selectedRun> --ellipse_multiplicative_factor <ellipse_multiplicative_factor> --precisionFactor <precisionFactor>"
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "ellipse_multiplicative_factor=",
                                                                 "precisionFactor="])
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules_plots_pixelConversion.py %s'%str1
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules_plots_pixelConversion.py %s'%str1
            sys.exit()
        elif option == "--selectedRun":
            runNumber = value.zfill(4)
        elif option == "--ellipse_multiplicative_factor":
            ellipse_multiplicative_factor = float(value)
        elif option == "--precisionFactor":
            precision_factor = int(value)
    
    warnings.filterwarnings("ignore")
    
    inputFolder = './Output_r%s/Output_imageSums_moduleDisplacements'%runNumber
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/Figures'%runNumber
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
    print 'Plot superorbit (2, 11), summed sectors and Gauss fits with fixed sigmas.\n'
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
                    bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
                    
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
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
                print bgSubtracted_total_sum.shape
                
                ax.set_title('%s'%(alpha_label[idx]), fontsize=11, loc='left')
                ax.set_title('{(%d, %d)}'%(h, k), fontsize=11)
                ax.imshow(bgSubtracted_total_sum, origin='lower', interpolation='nearest', vmin=0, vmax=3)
                ax.axis('off')
        
    matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i.png'%(outputFolder), dpi = 4*96 )       
    matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i.pdf'%(outputFolder), dpi = 4*96 )      
    matplotlib.pyplot.close()  
    
    
    # IMAGE SUMS SUBPLOTS 
    print 'Plot superorbit (2, 11), reconverted pixels.\n'
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
                bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
                bgSubtracted_total_sum_converted = convert(bgSubtracted_total_sum, precision_factor)
             
                ax.set_title('%s'%(alpha_label[idx]), fontsize=11, loc='left')
                ax.set_title('{(%d, %d)}'%(h, k), fontsize=11)
                ax.imshow(bgSubtracted_total_sum_converted, origin='lower', interpolation='nearest', vmin=0, vmax=3)
                ax.axis('off')
        
    matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i_converted.png'%(outputFolder), dpi = 4*96 )       
    matplotlib.pyplot.savefig('%s/Figure_bgsub_rotatetd_orbit_2_11_labels_f_to_i_converted.pdf'%(outputFolder), dpi = 4*96 )      
    matplotlib.pyplot.close()  
    
     
    
    
    
    ###
    orbit_label = [2, 11]
    h = orbit_label[0]
    k = orbit_label[1]
    print 'Plot Gauss fit and integration area.\n'
    print 'Plotting orbit (%d, %d).'%(h, k)
    
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
            
            data_toFit = data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
            gaussFit   = data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
            
            # SUMMED IMAGE WITH GAUSS FIT AND ELLIPTICAL INTEGRATION AREA
            matplotlib.pyplot.figure()
            ax = matplotlib.pyplot.gca()        
            ax.imshow(data_toFit, origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                        numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                        gaussFit, 4, linestyles='dashed', colors='w', linewidth=3)
            ax.axis('off')
            matplotlib.pyplot.plot((refined_x0, refined_x0), (0, bgSubtracted_total_sum.shape[0]), color='red', linestyle='dashed', linewidth=2)
            matplotlib.pyplot.plot((0, bgSubtracted_total_sum.shape[1]), (refined_y0, refined_y0), color='red', linestyle='dashed', linewidth=2)
            
            mean    = [bgSubtracted_total_sum.shape[1]/2 ,  bgSubtracted_total_sum.shape[0]/2]
            width   = 2* (ellipse_multiplicative_factor*sigma_x*10**precision_factor)
            height  = 2* (ellipse_multiplicative_factor*sigma_y*10**precision_factor)
            angle   = 0
            ellipse = matplotlib.patches.Ellipse(xy=mean, width=width, height=height, angle=angle, ec='m', fc='none', linewidth=3)
            ax.add_patch(ellipse)
        
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
            matplotlib.pyplot.close()  
            print 'Gaussian center: ', refined_x0, refined_y0
            
            # SUMMED IMAGE WITH GAUSS FIT AND ELLIPTICAL INTEGRATION AREA, CONVERTED PIXELS
            print 'Plotting image with converted pixels.'
            matplotlib.pyplot.figure()
            ax = matplotlib.pyplot.gca()    
            data_toFit_converted = convert(data_toFit, precision_factor)
            ax.imshow(data_toFit_converted, origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                        numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                        gaussFit, 4, linestyles='dashed', colors='w', linewidth=3)
            ax.axis('off')
            matplotlib.pyplot.plot((refined_x0, refined_x0), (0, bgSubtracted_total_sum.shape[0]), color='red', linestyle='dashed', linewidth=2)
            matplotlib.pyplot.plot((0, bgSubtracted_total_sum.shape[1]), (refined_y0, refined_y0), color='red', linestyle='dashed', linewidth=2)
            
            mean    = [bgSubtracted_total_sum.shape[1]/2 ,  bgSubtracted_total_sum.shape[0]/2]
            width   = 2* (ellipse_multiplicative_factor*sigma_x*10**precision_factor)
            height  = 2* (ellipse_multiplicative_factor*sigma_y*10**precision_factor)
            angle   = 0
            ellipse = matplotlib.patches.Ellipse(xy=mean, width=width, height=height, angle=angle, ec='m', fc='none', linewidth=3)
            ax.add_patch(ellipse)
        
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d_converted.png'%(outputFolder, h, k), dpi = 4*96 )       
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
            matplotlib.pyplot.close()  
            print 'Gaussian center: ', refined_x0, refined_y0
            
            # HORIZONTAL EXTRACT (DATA AND FIT) 
            print 'Plotting horizontal extract.'
            horizontal_extract = data_toFit[int(refined_y0), :]
            intervals = numpy.linspace(0, len(horizontal_extract), 21)
            avgs = []
            ns = []
            for i in range(0, len(intervals)-1):
                n1 = int(intervals[i])
                n2 = int(intervals[i+1])
                ptsToAvg = horizontal_extract[n1:n2]
                avg = numpy.average(ptsToAvg)
                avgs.append(avg)
                n = (float(n1+n2))/2
                ns.append(n)
                
            matplotlib.pyplot.figure()
            matplotlib.pyplot.scatter(ns, avgs, edgecolor='none')
            horizontal_extract_fit = gaussFit[int(refined_y0), :]
            matplotlib.pyplot.plot(horizontal_extract_fit, 'c', linewidth=2)
            matplotlib.pyplot.gca().set_xlim([0, 200])
            matplotlib.pyplot.gca().set_xticklabels([])
            matplotlib.pyplot.gca().set_ylabel(r"$I$ (photons)", fontsize = 22, rotation = 'vertical')
            matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=10)   
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig('%s/Gauss_Fit_horizontalExtract_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
            matplotlib.pyplot.savefig('%s/Gauss_Fit_horizontalExtract_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
            matplotlib.pyplot.close()  
            
            # VERTICAL EXTRACT (DATA AND FIT)
            print 'Plotting vertical extract.'
            vertical_extract = data_toFit[:, int(refined_x0)]
            intervals = numpy.linspace(0, len(vertical_extract), 21)
            avgs = []
            ns = []
            for i in range(0, len(intervals)-1):
                n1 = int(intervals[i])
                n2 = int(intervals[i+1])
                ptsToAvg = vertical_extract[n1:n2]
                avg = numpy.average(ptsToAvg)
                avgs.append(avg)
                n = (float(n1+n2))/2
                ns.append(n)
                
            matplotlib.pyplot.figure()
            matplotlib.pyplot.scatter(ns, avgs, edgecolor='none')
            vertical_extract_fit = gaussFit[:, int(refined_x0)]
            matplotlib.pyplot.plot(vertical_extract_fit, 'c', linewidth=2)
            matplotlib.pyplot.gca().set_xlim([0, 200])
            matplotlib.pyplot.gca().set_xticklabels([])
            matplotlib.pyplot.gca().set_ylabel(r"$I$ (photons)", fontsize = 22, rotation = 'vertical')
            matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=10)    
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig('%s/Gauss_Fit_verticalExtract_%d_%d.png'%(outputFolder, h, k), dpi = 4*96 )       
            matplotlib.pyplot.savefig('%s/Gauss_Fit_verticalExtract_%d_%d.pdf'%(outputFolder, h, k), dpi = 4*96 )      
            matplotlib.pyplot.close()
    
    flag = 0
    if flag ==1:
        print 'Plotting all orbits.'
        for key, orbit in orbitsDictionary.items():
                
            label = orbit.label
            h_label = label[0]
            k_label = label[1]
                    
            bgSubtracted_total_sum = orbit.bgSubtracted_total_sum[10**precision_factor*5 : 10**precision_factor*25, 10**precision_factor*5 : 10**precision_factor*25]
            
            q = 2 * numpy.pi / orbit.resolution        
            sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
            sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2]) 
            
            refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit_fixed_sigmas(bgSubtracted_total_sum, 
                                                                                                                                sigma_x*10**precision_factor, 
                                                                                                                                sigma_y*10**precision_factor)
            
            data_toFit = data.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
            gaussFit   = data_fitted.reshape(bgSubtracted_total_sum.shape[0], bgSubtracted_total_sum.shape[1])
            
            # SUMMED IMAGE WITH GAUSS FIT AND ELLIPTICAL INTEGRATION AREA
            matplotlib.pyplot.figure()
            ax = matplotlib.pyplot.gca()        
            ax.imshow(data_toFit, origin='lower', interpolation='nearest', vmin=0, vmax=3)
            ax.contour(numpy.linspace(0, bgSubtracted_total_sum.shape[0]-1, bgSubtracted_total_sum.shape[0]), 
                                        numpy.linspace(0, bgSubtracted_total_sum.shape[1]-1, bgSubtracted_total_sum.shape[1]), 
                                        gaussFit, 4, linestyles='dashed', colors='w', linewidth=3)
            ax.axis('off')
            matplotlib.pyplot.plot((refined_x0, refined_x0), (0, bgSubtracted_total_sum.shape[0]), color='red', linestyle='dashed', linewidth=2)
            matplotlib.pyplot.plot((0, bgSubtracted_total_sum.shape[1]), (refined_y0, refined_y0), color='red', linestyle='dashed', linewidth=2)
            
            mean    = [bgSubtracted_total_sum.shape[1]/2 ,  bgSubtracted_total_sum.shape[0]/2]
            width   = 2* (ellipse_multiplicative_factor*sigma_x*10**precision_factor)
            height  = 2* (ellipse_multiplicative_factor*sigma_y*10**precision_factor)
            angle   = 0
            ellipse = matplotlib.patches.Ellipse(xy=mean, width=width, height=height, angle=angle, ec='m', fc='none', linewidth=3)
            ax.add_patch(ellipse)
        
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.png'%(outputFolder, h_label, k_label), dpi = 4*96 )       
            matplotlib.pyplot.savefig('%s/Gauss_Fit_%d_%d.pdf'%(outputFolder, h_label, k_label), dpi = 4*96 )      
            matplotlib.pyplot.close()  
    print 'Done.'

if __name__ == "__main__":
    print "\n**** CALLING imageSums_displaced_modules_plots_pixelConversion ****"
    plot(sys.argv[1:])