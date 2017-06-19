# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import makeOrbits
import imageSums_utilities
import detectorModules

import buildReciprocalLattice
import calculatePredictedPattern


            
def calculate_moduleDisplacements_Function():
    
    # PARAMETERS
    selectedRun = '0127'
    resolutionLimit = 4.0
    intensity_threshold = 2
    nCountsPerPhoton = 26
    halfWidth = 25

    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements'%selectedRun    
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
            
    # PREPARE LATTICES TO IMAGES MATCHING
    imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%selectedRun    
    lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                           %(selectedRun, selectedRun, selectedRun)) # 629 lattices  
    lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(selectedRun, selectedRun, selectedRun), 'r')
    lattices_names = list(lattices_names)    
    images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(selectedRun, selectedRun))
    images_names = list(images_names)
    
    # EXTRACT GEOMETRY
    geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5'     # same for all runs
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    
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
    
    print xGeometry_np.shape, yGeometry_np.shape
    print row_interface
    print column_interface

    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                    # 220 orbits to 4 A
    print '%d orbits'%len(orbits)
        
    # CALCULATE STANDARD PATTERN    
    nominalCell = 62.45
    hmax = 100
    kmax = 100
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(nominalCell, hmax, kmax, resolutionLimit)
    tiltAngle = 0
    trialInPlaneRotation = 0
    wavelength = 1.485 # A
    wavevector = 2 * numpy.pi/wavelength
    detectorDistance = 0.235 # m
    pixelSize = 0.000110 # m
    standardPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, wavevector, tiltAngle, detectorDistance, pixelSize)
        
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
            print '\n*** MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound
            
            # CALCULATE MODULE RESOLUTION RANGE (speed up X4)
            border = 5
            bottom = bottom_bound + border
            top = top_bound - border
            left = left_bound + border
            right = right_bound - border
            
            x_left_bottom = xGeometry_np[bottom, left]   
            y_left_bottom = yGeometry_np[bottom, left]
            d_left_bottom = numpy.sqrt(x_left_bottom**2 + y_left_bottom**2)    # m
            
            x_right_bottom = xGeometry_np[bottom, right]   
            y_right_bottom = yGeometry_np[bottom, right]
            d_right_bottom = numpy.sqrt(x_right_bottom**2 + y_right_bottom**2) # m
            
            x_left_top = xGeometry_np[top, left]   
            y_left_top = yGeometry_np[top, left]
            d_left_top = numpy.sqrt(x_left_top**2 + y_left_top**2)             # m
            
            x_right_top = xGeometry_np[top, right]   
            y_right_top = yGeometry_np[top, right]
            d_right_top = numpy.sqrt(x_right_top**2 + y_right_top**2)          # m
            
            # RELAXED DISTANCES
            minDistance = min([d_left_bottom, d_left_top, d_right_bottom, d_right_top]) - 0.01   # m
            maxDistance = max([d_left_bottom, d_left_top, d_right_bottom, d_right_top]) + 0.01   # m
            print 'MODULE DISTANCE FROM DIRECT BEAM RANGE: [%f, %f]'%(minDistance, maxDistance)
            
            moduleFolder = '%s/Module_%d_%d'%(outputFolder, top_bound, right_bound)
            if not os.path.exists(moduleFolder):
                os.mkdir(moduleFolder)
                
            module_x0s = []
            module_y0s = []
            
            # LOOP ON ORBITS
            for orbit in orbits:
        
                indices_1 = orbit.orbitIndices[0]
                indices_2 = orbit.orbitIndices[1]
                indices_3 = orbit.orbitIndices[2]
                h_1 = indices_1[0]
                k_1 = indices_1[1]
                h_2 = indices_2[0]
                k_2 = indices_2[1]
                h_3 = indices_3[0]
                k_3 = indices_3[1]
                
                # EXTRACT ORBIT LABEL:
                label = orbit.label
                h_label = label[0]
                k_label = label[1]  
                
                print '\nOrbit %d, %d'%(h_label, k_label)       
                
                # CHECK ORBIT INTENSITY AND DETECTOR DISTANCE
                orbit_use_flag = 0
                orbit_use_flag_I = 0
                orbit_use_flag_D = 0
                
                I = numpy.nan
                spotDistance = numpy.nan
                
                for rough_intensity in rough_intensities_table:
                    if h_label == rough_intensity[0] and k_label == rough_intensity[1]:
                        I = float(rough_intensity[2])
                        if I >= intensity_threshold*nCountsPerPhoton:
                            orbit_use_flag_I = 1
                            
                for predictedSpot in standardPattern:
                    if predictedSpot[0] == h_label and predictedSpot[1] == k_label:
                        spotDistance = float(predictedSpot[10]) * pixelSize # m
                        if minDistance <= spotDistance <= maxDistance:
                            orbit_use_flag_D = 1
                            
                if orbit_use_flag_I == 1 and orbit_use_flag_D == 1:
                    orbit_use_flag = 1
                    
                ###
                if orbit_use_flag == 0:
                    print 'NOT using orbit %d, %d'%(h_label, k_label) 
                    print 'I = ', I
                    print minDistance, spotDistance, maxDistance
                ###
                    
                else:
                    print 'Using orbit %d, %d'%(h_label, k_label)
                    print 'I = ', I
                    print minDistance, spotDistance, maxDistance
                    
                    # USE THIS ORBIT TO CALCULATE MODULE DISPLACEMENT
                    module_sum = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 50x50
                    module_sum = numpy.matrix(module_sum)    
                    nTerms_module_sum = 0
                   
                    # LOOP ON ALL PROCESSED LATTICES (FIXED ORBIT, FIXED MODULE)
                    for index in range(0, len(lattices_names)):
                        lattice_name = lattices_names[index]
                        imageNumber = int(lattice_name[80:84])
                        latticeMatrix = lattices[index]                    # h_t k_t q_rod I_scaled flag i_unassembled j_unassembled scale
                        image_name = images_names[imageNumber-1].split()[1]
                        
                        if latticeMatrix[0, 4] == 0:                       # Check transformation and scaling flag
                            continue
                       
                        # LOAD UNASSEMBLED IMAGE
                        unassembledDataFile = h5py.File('%s/%s'%(imagesDirectoryName, image_name), 'r')
                        unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
                        unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####  (1480, 1552)
                        nRows = unassembledData.shape[0]
                        nColumns = unassembledData.shape[1]
                        
                        # IN THE CURRENT LATTICE, LOOP ON SPOTS BELONGING TO THE CURRENT ORBIT  
                        for spot in latticeMatrix:
                            h = int(spot[0])
                            k = int(spot[1])
                            if (h == h_1 and k == k_1) or (h == h_2 and k == k_2) or (h == h_3 and k == k_3):
                                i = int(spot[5])
                                j = int(spot[6])
                                
                                # CHECK MODULE
                                if not (bottom_bound <= i < top_bound):
                                    continue
                                if not (left_bound <= j < right_bound):
                                    continue
                                
                                left_edge   = j - halfWidth
                                right_edge  = j + halfWidth
                                bottom_edge = i - halfWidth
                                top_edge    = i + halfWidth
                                if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                                    continue
                                
                                ### EXTRACT 50x50 SECTOR ###
                                spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]  # 50x50
                                
                                ### SCALE ###
                                spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                                spotMatrix = spot[7] * spotMatrix
                                
                                ### DETERMINE AZIMUTH ON DETECTOR ###
                                detectorAzimuth = imageSums_utilities.calculate_detectorAzimuth(xGeometry_np, yGeometry_np, i, j)
                                        
                                ### DETERMINE MODULE ROTATION ANGLE ###
                                moduleRotation = imageSums_utilities.calculate_moduleRotation(xGeometry_np, yGeometry_np, i, j)
                                 
                                ### ROTATION ANGLE ###
                                rotationAngle = - detectorAzimuth + moduleRotation
                                rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                                
                                ### ROTATE ###  ### CLOCKWISE !!! ###
                                spotMatrix_rotated = imageSums_utilities.clockWiseRotation(spotMatrix, rotationAngle)
                                
                                ### SUM ###
                                module_sum = module_sum + spotMatrix_rotated # SUM ON SINGLE ORBIT AND SINGLE MODULE                                
                                nTerms_module_sum = nTerms_module_sum + 1
                                                    
                    print 'nTerms_module: ', nTerms_module_sum
                    if nTerms_module_sum == 0:
                        continue
                                                                
                    ### CALCULATE BACKGROUND ###
                    background = imageSums_utilities.calculateBackground_noImg(module_sum)
                              
                    ### BG SUBTRACTION ###
                    bgSubtracted_module_sum = module_sum - background
                    
                    ### NORMALIZATION ###
                    bgSubtracted_module_sum_normalized = bgSubtracted_module_sum / nTerms_module_sum
                    
                    ### SINGLE MODULE PLOT ###
                    myFigureObject = matplotlib.pyplot.figure()
                    myAxesImageObject = matplotlib.pyplot.imshow(bgSubtracted_module_sum_normalized, origin='lower', interpolation='nearest')
                    matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f'%(h_label, k_label, orbit.resolution))
                    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)                    
                    matplotlib.pyplot.savefig('%s/bgsub_rotatetd_top_%d_right_%d_n_%d_orbit_%d_%d.png'
                                               %(moduleFolder, top_bound, right_bound, nTerms_module_sum, h_label, k_label), dpi = 2*96 )                    
                    matplotlib.pyplot.close()  
                    
                    ### TRY GAUSSIAN FIT OF BG SUBTRACTED, ROTATED SUM ON SINGLE MODULE ###
                    refined_sigma_x, refined_sigma_y, refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted = imageSums_utilities.do_gaussFit(bgSubtracted_module_sum_normalized)
                    
                    if numpy.isnan(gauss_integral):
                        print "Gauss integral isnan"
                    else:
                        # SELECT GOOD FIT         
                        if nTerms_module_sum >= 50 and abs(gauss_integral) > intensity_threshold*nCountsPerPhoton and refined_amplitude > 0 and \
                        refined_sigma_x < 4 and refined_sigma_y < 4:
                            
                            module_x0s.append(refined_x0)
                            module_y0s.append(refined_y0)
                            
                            ### PLOT GAUSS FIT, SINGLE MODULE ###
                            myFigureObject = matplotlib.pyplot.figure()
                            myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*halfWidth, 2*halfWidth), origin='lower', interpolation='nearest')
                            matplotlib.pyplot.gca().contour(numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                            numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                            data_fitted.reshape(2*halfWidth, 2*halfWidth), 4, colors='w')
                            matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f Gauss integral: %.1f counts'
                                                     %(h_label, k_label, orbit.resolution, gauss_integral))
                            myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                            matplotlib.pyplot.savefig('%s/gauss_bgsub_rotatetd_top_%d_right_%d_n_%d_orbit_%d_%d.png'
                                                       %(moduleFolder, top_bound, right_bound, nTerms_module_sum, h_label, k_label), dpi = 2*96 )       
                            matplotlib.pyplot.close()  

            print module_x0s
            print module_y0s
            avg_x0 = numpy.average(module_x0s)
            avg_y0 = numpy.average(module_y0s)
            std_x0 = numpy.std(module_x0s)
            std_y0 = numpy.std(module_y0s)
            print avg_x0, std_x0, avg_y0, std_y0
            x0s_file = open('%s/module_%d_%d_x0s.pkl'%(outputFolder, top_bound, right_bound), 'wb')
            pickle.dump(module_x0s, x0s_file)
            x0s_file.close()
            y0s_file = open('%s/module_%d_%d_y0s.pkl'%(outputFolder, top_bound, right_bound), 'wb')
            pickle.dump(module_y0s, y0s_file)
            y0s_file.close()

if __name__ == "__main__":
    print "\n**** CALLING calculate_moduleDisplacements ****"
    calculate_moduleDisplacements_Function()