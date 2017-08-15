# -*- coding: utf-8 -*-
import numpy
import os
import pickle
import sys
import getopt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import makeOrbits
import imageSums_utilities


            
def imageSums(myArguments):

    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "resolutionLimit=", 
                                                                 "halfWidth=",
                                                                 "lowFluctuationThreshold=",
                                                                 "precisionFactor="])
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules.py --selectedRun <selectedRun> --resolutionLimit <resolutionLimit> --halfWidth <halfWidth> --lowFluctuationThreshold <lowFluctuationThreshold> --precisionFactor <precisionFactor>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules.py --selectedRun <selectedRun> --resolutionLimit <resolutionLimit> --halfWidth <halfWidth> --lowFluctuationThreshold <lowFluctuationThreshold> --precisionFactor <precisionFactor>'
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--lowFluctuationThreshold":
            bg_lowFluctuationThreshold = float(value)
        elif option == "--precisionFactor":
            precisionFactor = int(value)
            
    truncated_halfWidth = halfWidth - 10
        
    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements'%selectedRun    
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
    if not os.path.exists('%s/high_res_7_6'%outputFolder):
        os.mkdir('%s/high_res_7_6'%outputFolder)
    if not os.path.exists('%s/high_res_6_5'%outputFolder):
        os.mkdir('%s/high_res_6_5'%outputFolder)
    if not os.path.exists('%s/high_res_5_4'%outputFolder):
        os.mkdir('%s/high_res_5_4'%outputFolder)
    if not os.path.exists('%s/low_res'%outputFolder):
        os.mkdir('%s/low_res'%outputFolder)
            
    # LOGGING
    fOpen = open('%s/imageSums_displaced_modules.txt'%outputFolder, 'w')
    fOpen.write('h k Gauss_amplitude x0 y0 sigma_x sigma_y I_Gauss I_sum\n')
   
    # EXTRACT MODULE DISPLACEMENTS
    module_displacements_file = open('./Output_r%s/ModuleDisplacements/module_displacements.pkl'%selectedRun, 'rb')
    modules = pickle.load(module_displacements_file) # bottom, top, left, right, N, <Dx0>, sig_Dx0, <Dy0>, sig_Dy0
    module_displacements_file.close()
    print 'MODULES SHAPE: ', modules.shape
    
    # EXTRACT PARTIAL SUMS
    partialSums_file = open('./Output_r%s/ModuleDisplacements/partialSums_list.pkl'%selectedRun, 'rb')
    partialSums_list = pickle.load(partialSums_file)
    partialSums_file.close()
    print '%d partial sums'%len(partialSums_list)
    
    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                   # 220 orbits to 4 A
    print '%d orbits'%len(orbits)
       
    # FOR EACH ORBIT, DO SUMS OF SINGLE-MODULE SUMS
    orbits_withSum = {}
    for orbit in orbits:
        
        # EXTRACT ORBIT LABEL:
        label = orbit.label
        orbit_h_label = label[0]
        orbit_k_label = label[1]                           
        print 'ORBIT: ', orbit.label 
        
        # FIGURES FOLDER
        if orbit.resolution > 7.0:
            orbit_figureFolder = '%s/low_res'%outputFolder                                      
        if 6.0 < orbit.resolution <= 7.0:
            orbit_figureFolder = '%s/high_res_7_6'%outputFolder            
        if 5.0 < orbit.resolution <= 6.0:
            orbit_figureFolder = '%s/high_res_6_5'%outputFolder            
        if 4.0 < orbit.resolution <= 5.0:
            orbit_figureFolder = '%s/high_res_5_4'%outputFolder
                     
        # ONE ORBIT
        nTot_orbit = 0
        total_sum = numpy.zeros(( 2*truncated_halfWidth*(10**precisionFactor), 2*truncated_halfWidth*(10**precisionFactor) ))
        
        # LOOP ON PARTIAL SUMS
        for i in range(0, len(partialSums_list)):                       
            partialSum = partialSums_list[i]
            if (partialSum.h_label == orbit_h_label and 
                partialSum.k_label == orbit_k_label):
                                        
                partialSum_sector = partialSum.partialSum # NOT BG SUBTRACTED, NOT NORMALIZED ON N OF TERMS, EXPRESSED IN COUNTS
                nTerms = partialSum.nTerms   
                
                bottomBound = partialSum.module_bottomBound
                topBound    = partialSum.module_topBound
                leftBound   = partialSum.module_leftBound 
                rightBound  = partialSum.module_rightBound 
                    
                module_Dx0 = numpy.nan
                module_Dy0 = numpy.nan
                
                for module_line in modules:
                    if (module_line[0] == bottomBound and
                        module_line[1] == topBound and
                        module_line[2] == leftBound and
                        module_line[3] == rightBound):
                        module_Dx0 = module_line[5]
                        module_Dy0 = module_line[7]
                
                if (numpy.isnan(module_Dx0) or numpy.isnan(module_Dy0)):
                    print 'PROBLEM'
                    
                module_x0 = halfWidth + module_Dx0
                module_y0 = halfWidth + module_Dy0
                    
                recentered_partialSum_sector = imageSums_utilities.recenter(partialSum_sector, 
                                                                             module_x0, 
                                                                             module_y0,
                                                                             precisionFactor,
                                                                             truncated_halfWidth) # 50x50 -> 300x300, NO BG SUB, NO NORMALIZATION, DETECTOR COUNTS
                total_sum = total_sum + recentered_partialSum_sector
                nTot_orbit = nTot_orbit + nTerms                                                        
                                                        
        # NORMALIZATION
        total_sum = total_sum / nTot_orbit 
        
        # BG SUBTRACTION
        bg, n_bgPixels = imageSums_utilities.calculateBackground_nBg(total_sum, lowFluctuationThreshold=bg_lowFluctuationThreshold)
        bgSubtracted_total_sum = total_sum - bg
        
        # GENERATE UPDATED ORBITS DICTIONARY
        orbit.nTerms = nTot_orbit
        orbit.n_bgPixels = n_bgPixels
        orbit.bgSubtracted_total_sum = bgSubtracted_total_sum
        orbits_withSum['%d_%d'%(orbit_h_label, orbit_k_label)] = orbit
        
        ### SINGLE PLOT ###
        myFigureObject = matplotlib.pyplot.figure()
        myAxesImageObject = matplotlib.pyplot.imshow(bgSubtracted_total_sum, origin='lower', interpolation='nearest')
        matplotlib.pyplot.title('Orbit: %d %d'%(orbit_h_label, orbit_k_label))
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)        
        matplotlib.pyplot.savefig('%s/modules_sum_%d_%d_nTot_%d.png'
                                   %(orbit_figureFolder, orbit_h_label, orbit_k_label, nTot_orbit), dpi = 2*96 )
        matplotlib.pyplot.close()  
        
        ### SUM INTEGRATION ###
        integratedIntensity = imageSums_utilities.integrate(bgSubtracted_total_sum, expansion_factor=precisionFactor)  
        integratedIntensity = integratedIntensity/((10**precisionFactor)**2)
        
        ### GAUSS FIT ###
        refined_sigma_x, refined_sigma_y, \
        refined_x0, refined_y0, \
        refined_amplitude, gauss_integral, \
        data, data_fitted = imageSums_utilities.do_gaussFit(bgSubtracted_total_sum)
        
        if not numpy.isnan(gauss_integral):
            refined_sigma_x = refined_sigma_x/(10**precisionFactor) # Can be negative at this stage.
            refined_sigma_y = refined_sigma_y/(10**precisionFactor)
            gauss_integral  = gauss_integral/((10**precisionFactor)**2)
            
        print integratedIntensity, gauss_integral
        
        ### PLOT GAUSS FIT ###
        if not numpy.isnan(gauss_integral):
            
            myFigureObject = matplotlib.pyplot.figure()
            myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*truncated_halfWidth*(10**precisionFactor), 2*truncated_halfWidth*(10**precisionFactor)), 
                                                         origin='lower', interpolation='nearest')
            try:
                matplotlib.pyplot.gca().contour(numpy.linspace(0, 2*truncated_halfWidth*(10**precisionFactor)-1, 2*truncated_halfWidth*(10**precisionFactor)), 
                                                numpy.linspace(0, 2*truncated_halfWidth*(10**precisionFactor)-1, 2*truncated_halfWidth*(10**precisionFactor)), 
                                                data_fitted.reshape(2*truncated_halfWidth*(10**precisionFactor), 2*truncated_halfWidth*(10**precisionFactor)), 4, colors='w')
   
            except:
                print 'PROBLEM DRAWING CONTOURS'
                
            matplotlib.pyplot.title('Orbit: %d %d Gauss integral: %.1f counts\nSig_x %.2f Sig_y %.2f'
                                     %(orbit_h_label, orbit_k_label, gauss_integral, refined_sigma_x, refined_sigma_y))
            myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
            matplotlib.pyplot.savefig('%s/gauss_module_sum_%d_%d_nTot_%d.png'
                                       %(orbit_figureFolder, orbit_h_label, orbit_k_label, nTot_orbit), dpi = 2*96)
            matplotlib.pyplot.close()  
            
        fOpen.write('%4d%4d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n'%(orbit_h_label, orbit_k_label, 
                                                                   refined_amplitude, refined_x0, refined_y0,
                                                                   refined_sigma_x, refined_sigma_y,
                                                                   gauss_integral, integratedIntensity))
    fOpen.close()
    
    # SAVE ORBIT OBJECTS
    allOrbits_file = open('%s/orbits.pkl'%outputFolder, 'wb')
    pickle.dump(orbits_withSum, allOrbits_file)
    allOrbits_file.close()    
    


if __name__ == "__main__":
    print "\n**** CALLING imageSums_displaced_modules ****"
    imageSums(sys.argv[1:])   