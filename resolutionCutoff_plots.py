# -*- coding: utf-8 -*-
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle

import imageSums_utilities

def resolutionCutoff_plots():
    
    # PARAMETERS
    nCountsPerPhoton = 26
    multiplicative_factor = 2.5
    
    # FOLDERS
    outputFolder = './Output_resolution_cutoff'
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)       
        
### CALCULATE I/SIG(I) IN THE IMAGE SUMS METHOD

    # EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
    fRead = open('./Output_imageSums_sigmaFits/sigmaXCurveParameters.pkl', 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()
    
    fRead = open('./Output_imageSums_sigmaFits/sigmaYCurveParameters.pkl', 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close()   
    
    # LOAD ORBITS
    fileToOpen = './Output_imageSums/orbits.pkl'
    fRead = open(fileToOpen, 'rb')
    orbitsDictionary = pickle.load(fRead)                                               
    fRead.close()
    
    # CHECK HOW MANY ORBITS ARE THERE (220)
    nOrbits = 0
    for key, orbit in orbitsDictionary.items():
        nOrbits = nOrbits + 1    
    print nOrbits
    
    # LOOP ON ORBITS
    qs = []
    ratios = []
    for key, orbit in orbitsDictionary.items():
        
        label = orbit.label
        h_label = label[0]
        k_label = label[1]
        print h_label, k_label
        bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtractedOrbitSumMatrix_rotated
        q = 2 * numpy.pi / orbit.resolution
        sigma_x = imageSums_utilities.quadratic(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.line_plus_sigmoid(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2], sigmaYCurveParameters[3], sigmaYCurveParameters[4])    
            
        # N COUNTS TO N PHOTONS CONVERSION
        bgSubtractedOrbitSumMatrix_rotated =  bgSubtractedOrbitSumMatrix_rotated/nCountsPerPhoton
        
        # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
        integrationMask = numpy.zeros((bgSubtractedOrbitSumMatrix_rotated.shape))   
        ringMask = numpy.zeros((bgSubtractedOrbitSumMatrix_rotated.shape)) 
        
        colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
        x_axis = multiplicative_factor * sigma_x
        y_axis = multiplicative_factor * sigma_y
        centerX = integrationMask.shape[1] / 2
        centerY = integrationMask.shape[0] / 2
        distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
        
        integrationMask[numpy.where(distance < 1)] = 1
        ringMask[numpy.where(distance > 1.5)] = 1
        ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
        
        # CONTROL PLOTS
        matplotlib.pyplot.imshow(integrationMask)
        matplotlib.pyplot.savefig('%s/integrationMask_%d_%d'%(outputFolder, h_label, k_label))
        matplotlib.pyplot.close()
        matplotlib.pyplot.imshow(ringMask)
        matplotlib.pyplot.savefig('%s/ringMask_%d_%d'%(outputFolder, h_label, k_label))
        matplotlib.pyplot.close()
        
        # CALCULATE RATIO
        n_peak = integrationMask.sum()
        n_bg = ringMask.sum()
        
        integratedIntensity = numpy.multiply(integrationMask, bgSubtractedOrbitSumMatrix_rotated).sum()
        avg_I = integratedIntensity/n_peak
           
        integratedBg = numpy.multiply(ringMask, bgSubtractedOrbitSumMatrix_rotated).sum()
        avg_bg = integratedBg/n_bg
        ring_bg = numpy.multiply(ringMask, bgSubtractedOrbitSumMatrix_rotated)
        ring_bg = ring_bg.flatten().T
        
        sum_sq = 0
        n = 0
        for bg_pxl in ring_bg:
            if bg_pxl != 0:
                n = n+1
                diff = bg_pxl - avg_bg
                sq = diff**2
                sum_sq = sum_sq + sq
        if n_bg != n:
            print 'n_bg different from n'
        bg_fluctuations = numpy.sqrt(sum_sq/n_bg)
        ratio = avg_I/bg_fluctuations
    
        qs.append(q)
        ratios.append(ratio)
    
    # BINNING    
    qs_plot_imageSum = []
    ratios_plot_imageSum = [] 
    bins = numpy.linspace(min(qs), max(qs), 25)
    for i in range(0, len(bins)-1):
        left_q = bins[i]
        right_q = bins[i+1]
        q_bin     = [qs[i]     for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
        ratio_bin = [ratios[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
        print len(ratio_bin)
        q_avg_bin = numpy.average(q_bin)
        ratio_avg_bin = numpy.average(ratio_bin)
        qs_plot_imageSum.append(q_avg_bin)
        ratios_plot_imageSum.append(ratio_avg_bin)
                
### EXTRACT DATA FROM SINGLE IMAGE METHOD        
    qs_file = open('./Output_resolutionCutoff_singleImage/qs.pkl', 'rb')
    qs = pickle.load(qs_file)
    qs_file.close()    
    
    ratios_file = open('./Output_resolutionCutoff_singleImage/ratios.pkl', 'rb')
    ratios = pickle.load(ratios_file)
    ratios_file.close()    
    print len(qs), len(ratios)
    
    # BINNING
    qs_plot_single_image = []
    ratios_plot_single_image = [] 
    bins = numpy.linspace(min(qs), max(qs), 25)
    for i in range(0, len(bins)-1):
        left_q = bins[i]
        right_q = bins[i+1]
        q_bin     = [qs[i]     for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
        ratio_bin = [ratios[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
        q_avg_bin = numpy.average(q_bin)
        ratio_avg_bin = numpy.average(ratio_bin)
        qs_plot_single_image.append(q_avg_bin)
        ratios_plot_single_image.append(ratio_avg_bin)
    
### PLOT I/sig(I) FROM IMAGE SUM METHOD AND SINGLE IMAGE METHOD    
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(qs_plot_imageSum, ratios_plot_imageSum, s=50, facecolors='m', edgecolors='m', alpha = 1, label='Image sums')
    matplotlib.pyplot.scatter(qs_plot_single_image, ratios_plot_single_image, s=50, marker='x', c='b', label='Single images')
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$\frac{\rm I}{\sigma({\rm I})}$", fontsize = 12, rotation = 'vertical')
    matplotlib.pyplot.gca().set_ylim([0.05, 150])
    matplotlib.pyplot.gca().axhline(y=1, color='k')
    matplotlib.pyplot.gca().set_yscale('log')
    legend = matplotlib.pyplot.gca().legend(scatterpoints=1, fontsize=10)
    legend.draw_frame(False)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/IoverSigI_comparison.png'%outputFolder)



if __name__ == "__main__":
    print "\n**** CALLING resolutionCutoff_plots ****"
    resolutionCutoff_plots()