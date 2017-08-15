# -*- coding: utf-8 -*-
import h5py
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
import detectorModules


            
def calculate_moduleDisplacements_Function(myArguments):
    
    # DEFAULTS:
    selectedRun = '0127'
    resolutionLimit = 4.0
    intensity_threshold = 2
    nCountsPerPhoton = 26
    
    str1 = '--selectedRun <selectedRun> --resolutionLimit <resolutionLimit>'
    str2 = '--intensity_threshold <intensity_threshold> --nTerms_threshold <nTerms_threshold> --sigma_threshold <sigma_threshold> --distance_threshold <distance_threshold>'
    str3 = '--nCountsPerPhoton <nCountsPerPhoton> --geometryFile <geometryFile>'

    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "resolutionLimit=", 
                                                                 "intensity_threshold=", 
                                                                 "nTerms_threshold=",
                                                                 "sigma_threshold=",
                                                                 "distance_threshold=",
                                                                 "nCountsPerPhoton=", 
                                                                 "geometryFile="])
    except getopt.GetoptError:
        print 'Error Usage: python calculate_moduleDisplacements.py %s %s %s'%(str1, str2, str3)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_moduleDisplacements.py %s %s %s'%(str1, str2, str3)
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--intensity_threshold":
            intensity_threshold = float(value)
        elif option == "--nTerms_threshold":
            nTerms_threshold = int(value)
        elif option == "--sigma_threshold":
            sigma_threshold = float(value)
        elif option == "--distance_threshold":
            distance_threshold = float(value)
        elif option == "--nCountsPerPhoton":
            nCountsPerPhoton = float(value)
        elif option == "--geometryFile":
            geometryFile = value
            
    # FOLDERS
    outputFolder = './Output_r%s/ModuleDisplacements'%selectedRun    
    
    # LOG
    logFile = open('%s/calculate_moduleDisplacements.log'%outputFolder, 'w')
    logFile.write('module_top module_right h k x0 y0 sigmaX sigmaY gaussI gaussA nTerms\n')
            
    # EXTRACT GEOMETRY
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    print 'xGeometry SHAPE: ', xGeometry_np.shape
    print 'yGeometry SHAPE: ', yGeometry_np.shape
    
    # EXTRACT DETECTOR MODULES BOUNDARIES
    modules = numpy.zeros((xGeometry_np.shape[0], xGeometry_np.shape[1]), dtype=numpy.int)
    modules = detectorModules.discontinuityMask(xGeometry_np, yGeometry_np, modules)
    
    row_interface = []
    row_interface.append(0)
    for i in range(0, modules.shape[0]):
        if modules[i, 0] == 1:
            row_interface.append(i)
    row_interface.append(modules.shape[0]-1)     #[0, 184, 369, 554, 739, 924, 1109, 1294, 1479]
    
    column_interface = []
    column_interface.append(0)
    for j in range(0, modules.shape[1]):
        if modules[0, j] == 1:
            column_interface.append(j)
    column_interface.append(modules.shape[1]-1)  #[0, 193, 387, 581, 775, 969, 1163, 1357, 1551]
    
    print 'MODULE ROW INTERFACES: ', row_interface
    print 'MODULE COLUMN INTERFACES: ', column_interface

    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                    # 220 orbits to 4 A
    print '%d orbits'%len(orbits)
    
    # EXTRACT PARTIAL SUMS
    partialSums_file = open('%s/partialSums_list.pkl'%(outputFolder), 'rb')
    partialSums_list = pickle.load(partialSums_file)
    partialSums_file.close()
    print '%d partial sums'%len(partialSums_list)
    
    # EXTRACT APPROXIMATE ORBIT INTENSITIES
    rough_intensities = open('./Output_r%s/Output_imageSums/h_k_Isum_Igauss_sigX_sigY.txt'%selectedRun, 'r')
    rough_intensities_table = []
    for rough_intensity in rough_intensities:
        line = [int(rough_intensity.split()[0]), int(rough_intensity.split()[1]), float(rough_intensity.split()[3])] # h, k, I_gauss
        rough_intensities_table.append(line)
    rough_intensities_table = numpy.asarray(rough_intensities_table)
    rough_intensities.close()
        
    # LOOP ON MODULES
    for module_row in range(0, len(row_interface)-1):
        for module_column in range(0, len(column_interface)-1):
            bottom_bound = row_interface[module_row]
            top_bound = row_interface[module_row+1]
            left_bound = column_interface[module_column]
            right_bound = column_interface[module_column+1]
            
            print '\n*** MODULE: bottom=%8d top=%8d left=%8d right=%8d'%(bottom_bound, top_bound, left_bound, right_bound)
            
            moduleFolder = '%s/Module_%d_%d'%(outputFolder, top_bound, right_bound)
            if not os.path.exists(moduleFolder):
                os.mkdir(moduleFolder)
                
            module_x0s = []
            module_y0s = []
            
            # LOOP ON ORBITS
            for orbit in orbits:
                
                # EXTRACT ORBIT LABEL:
                label = orbit.label
                orbit_h_label = label[0]
                orbit_k_label = label[1]  
                
                print label
                
                orbitFlag = 0
                for rough_intensity in rough_intensities_table:
                    if orbit_h_label == rough_intensity[0] and orbit_k_label == rough_intensity[1]:
                        I = float(rough_intensity[2])
                        if I >= intensity_threshold*nCountsPerPhoton:
                            # USE THIS ORBIT 
                            orbitFlag = 1
                            
                if orbitFlag == 1:
                    # LOOP ON PARTIAL SUMS
                    for i in range(0, len(partialSums_list)):                       
                        partialSum = partialSums_list[i]
                        if (partialSum.h_label == orbit_h_label and 
                            partialSum.k_label == orbit_k_label and
                            partialSum.module_bottomBound == bottom_bound and
                            partialSum.module_topBound == top_bound and
                            partialSum.module_leftBound == left_bound and
                            partialSum.module_rightBound == right_bound):
                                
                            nTerms = partialSum.nTerms                                        
                            if nTerms >= nTerms_threshold:
                                partialSum_sector = partialSum.partialSum # NOT BG SUBTRACTED, NOT NORMALIZED ON N OF TERMS, EXPRESSED IN COUNTS
                                halfWidth = int(partialSum_sector.shape[0]/2)
                                
                                ### CALCULATE BACKGROUND ###
                                background = imageSums_utilities.calculateBackground_noImg(partialSum_sector)
                                          
                                ### BG SUBTRACTION ###
                                bgSubtracted_partial_sum = partialSum_sector - background
                                
                                ### NORMALIZATION ###
                                bgSubtracted_partial_sum_normalized = bgSubtracted_partial_sum / nTerms
                        
                                ### SINGLE MODULE PLOT ###
                                single_module_plot_flag = 0
                                if single_module_plot_flag == 1:
                                    myFigureObject = matplotlib.pyplot.figure()
                                    myAxesImageObject = matplotlib.pyplot.imshow(bgSubtracted_partial_sum_normalized, origin='lower', interpolation='nearest')
                                    matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f'%(orbit_h_label, orbit_k_label, orbit.resolution))
                                    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)                    
                                    matplotlib.pyplot.savefig('%s/partialSum_module_top_%d_right_%d_n_%d_orbit_%d_%d.png'
                                                               %(moduleFolder, top_bound, right_bound, nTerms, orbit_h_label, orbit_k_label), dpi = 2*96 )                    
                                    matplotlib.pyplot.close()  
                            
                           
                                ### TRY GAUSSIAN FIT OF BG SUBTRACTED, ROTATED SUM ON SINGLE MODULE ###
                                refined_sigma_x, refined_sigma_y, \
                                refined_x0, refined_y0, \
                                refined_amplitude, gauss_integral, \
                                data, data_fitted = imageSums_utilities.do_gaussFit(bgSubtracted_partial_sum_normalized)
    
                                logFile.write('%6d %6d %4d %4d %6.2f %6.2f %6.2f %6.2f %8.2f %8.2f %8d\n'
                                              %(top_bound, right_bound, 
                                                orbit_h_label, orbit_k_label, 
                                                refined_x0, refined_y0, 
                                                refined_sigma_x, refined_sigma_y,
                                                gauss_integral, refined_amplitude,
                                                nTerms))
                                
                                # SELECT GOOD FIT         
                                if abs(gauss_integral) > intensity_threshold*nCountsPerPhoton and \
                                   refined_amplitude > 0 and \
                                   abs(refined_sigma_x) < sigma_threshold and \
                                   abs(refined_sigma_y) < sigma_threshold and \
                                   abs(refined_x0 - halfWidth) < distance_threshold and \
                                   abs(refined_y0 - halfWidth) < distance_threshold:
                                    
                                    module_x0s.append(refined_x0)
                                    module_y0s.append(refined_y0)
                                    
                                    ### PLOT GAUSS FIT, SINGLE MODULE ###
                                    myFigureObject = matplotlib.pyplot.figure()
                                    myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*halfWidth, 2*halfWidth), origin='lower', interpolation='nearest')
                                    matplotlib.pyplot.gca().contour(numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                                    numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                                    data_fitted.reshape(2*halfWidth, 2*halfWidth), 4, colors='w')
                                    matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f Gauss integral: %.1f counts'
                                                             %(orbit_h_label, orbit_k_label, orbit.resolution, gauss_integral))
                                    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                                    matplotlib.pyplot.savefig('%s/gauss_partialSum_top_%d_right_%d_n_%d_orbit_%d_%d.png'
                                                               %(moduleFolder, top_bound, right_bound, nTerms, orbit_h_label, orbit_k_label), dpi = 2*96 )       
                                    matplotlib.pyplot.close()  

            print 'MODULE X0S: ', module_x0s
            print 'MODULE Y0S: ', module_y0s
            avg_x0 = numpy.average(module_x0s)
            avg_y0 = numpy.average(module_y0s)
            std_x0 = numpy.std(module_x0s)
            std_y0 = numpy.std(module_y0s)
            print 'X0 = %f +- %f'%(avg_x0, std_x0)
            print 'Y0 = %f +- %f'%(avg_y0, std_y0)
            x0s_file = open('%s/module_%d_%d_x0s.pkl'%(outputFolder, top_bound, right_bound), 'wb')
            pickle.dump(module_x0s, x0s_file)
            x0s_file.close()
            y0s_file = open('%s/module_%d_%d_y0s.pkl'%(outputFolder, top_bound, right_bound), 'wb')
            pickle.dump(module_y0s, y0s_file)
            y0s_file.close()            
    

if __name__ == "__main__":
    print "\n**** CALLING calculate_moduleDisplacements ****"
    calculate_moduleDisplacements_Function(sys.argv[1:])   