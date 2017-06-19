# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle
import random
import sys
import getopt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import makeOrbits
import imageSums_utilities
import buildReciprocalLattice
import calculatePredictedPattern

           
def imageSums_CChalf(myArguments):
    
    # PARAMETERS
    selectedRun = '0127'
    resolutionLimit = 4.0
    halfWidth = 25
    truncated_halfWidth = 15
    precision_factor = 1
    ellipse_multiplicative_factor = 2.5
        
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["N_lattices="])
    except getopt.GetoptError:
        print 'python imageSums_displaced_modules_CChalf.py --N_lattices <N_lattices>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':            
            print 'python imageSums_displaced_modules_CChalf.py --N_lattices <N_lattices>'
            sys.exit()
        elif option == "--N_lattices":
            N_lattices = int(value)
            
    print 'N_lattices = ', N_lattices
    
    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_CChalf'%selectedRun
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
    
    # EXTRACT MODULE DISPLACEMENTS
    module_displacements_file = open('./Output_r%s/Output_imageSums_moduleDisplacements/module_displacements.pkl'%selectedRun, 'rb')
    modules = pickle.load(module_displacements_file) # bottom, top, left, right, N, <x0>, sig_x0, <y0>, sig_y0
    module_displacements_file.close()
    
    # ADD MODULE LIMITS
    modules_with_d_range = []
    for module in modules:
        bottom_bound = int(module[0])
        top_bound = int(module[1])
        left_bound = int(module[2])
        right_bound = int(module[3])
        print '\n*** MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound
        
        # CALCULATE MODULE RESOLUTION RANGE 
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
        
        minDistance = min([d_left_bottom, d_left_top, d_right_bottom, d_right_top])  - 0.01  # m
        maxDistance = max([d_left_bottom, d_left_top, d_right_bottom, d_right_top])  + 0.01  # m
        print 'MODULE DISTANCE FROM DIRECT BEAM RANGE: [%f, %f]'%(minDistance, maxDistance)
        
        module_with_d_range = [module[0], module[1], module[2], module[3], module[4], module[5], module[6], module[7], module[8], minDistance, maxDistance]
        modules_with_d_range.append(module_with_d_range)
        
    
    # EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
    fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaXCurveParameters.pkl'%selectedRun, 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()
    
    fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaYCurveParameters.pkl'%selectedRun, 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close()  
    
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
    
    for sampling in range(0, 10): 
        print sampling
        # LOGGING
        fOpen = open('%s/imageSums_displaced_modules_%d_lattices_halves_%d.txt'%(outputFolder, N_lattices, sampling), 'w')
        fOpen.write('# h k I1_ellipse_fixed_sigmas I2_ellipse_fixed_sigmas N1 N2\n')
    
        # MAKE ORBIT OBJECTS LIST
        orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                   # 220 orbits to 4 A
        
        # EXTRACT HALF THE LATTICES          
        random_sample = random.sample(range(N_lattices), N_lattices/2)
       
        # FOR EACH ORBIT, FOR EACH MODULE, DO SECTOR SUMS
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
            
            spotDistance = numpy.nan
            for predictedSpot in standardPattern:
                if predictedSpot[0] == h_label and predictedSpot[1] == k_label:
                    spotDistance = float(predictedSpot[10]) * pixelSize # m
                                   
            print 'ORBIT: ', orbit.orbitIndices, 'LABEL: ', orbit.label, 'DISTANCE', spotDistance
            
            # ONE ORBIT
            total_sum_1 = numpy.zeros((2*truncated_halfWidth*(10**precision_factor), 2*truncated_halfWidth*(10**precision_factor))) # 300x300
            total_sum_2 = numpy.zeros((2*truncated_halfWidth*(10**precision_factor), 2*truncated_halfWidth*(10**precision_factor))) # 300x300
            
            nTot_1_orbit = 0
            nUsed_1_orbit = 0
            
            nTot_2_orbit = 0
            nUsed_2_orbit = 0
            
            # LOOP ON ALL MODULES
            for module in modules_with_d_range:
                
                min_distance = float(module[9])
                max_distance = float(module[10])
                
                bottom_bound = int(module[0])
                top_bound = int(module[1])
                left_bound = int(module[2])
                right_bound = int(module[3])
                x0 = float(module[5])
                y0 = float(module[7])
                if numpy.isnan(x0):
                    x0 = float(halfWidth)
                if numpy.isnan(y0):
                    y0 = float(halfWidth)
                print 'MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound 
                print '<x0> = %.2f, <y0> = %.2f'%(x0, y0)
                
                # SKIP MODULES WITH RESOLUTION RANGE NOT INCLUDING THE CURRENT ORBIT
                if spotDistance < min_distance or spotDistance > max_distance:
                    print 'UNUSED MODULE [%.2f, %.2f] (ORBIT AT %.2f)'%(min_distance, max_distance, spotDistance)
                    continue
                
                module_sum_1 = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 50x50
                module_sum_1 = numpy.matrix(module_sum_1)
    
                nTerms_1_module = 0
                
                module_sum_2 = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 50x50
                module_sum_2 = numpy.matrix(module_sum_2)
    
                nTerms_2_module = 0
                
                # LOOP ON ALL PROCESSED LATTICES (FIXED ORBIT, FIXED MODULE)
                for index in range(0, N_lattices):
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
                            if index in random_sample:
                                module_sum_1 = module_sum_1 + spotMatrix_rotated # SUM ON SINGLE MODULE                        
                                nTerms_1_module = nTerms_1_module + 1
                                nTot_1_orbit = nTot_1_orbit + 1
                            else:
                                module_sum_2 = module_sum_2 + spotMatrix_rotated # SUM ON SINGLE MODULE                        
                                nTerms_2_module = nTerms_2_module + 1
                                nTot_2_orbit = nTot_2_orbit + 1
                                
                                
                print 'nTerms_1_module: ', nTerms_1_module
                print 'nTerms_2_module: ', nTerms_1_module
                if nTerms_1_module > 0:
                    recentered_sum_1 = imageSums_utilities.recenter(module_sum_1, x0, y0, precision_factor, truncated_halfWidth) # 50x50 -> 300x300, NO BG SUB, NO NORMALIZATION
                    total_sum_1 = total_sum_1 + recentered_sum_1
                    nUsed_1_orbit = nUsed_1_orbit + nTerms_1_module       
                if nTerms_2_module > 0:
                    recentered_sum_2 = imageSums_utilities.recenter(module_sum_2, x0, y0, precision_factor, truncated_halfWidth) # 50x50 -> 300x300, NO BG SUB, NO NORMALIZATION
                    total_sum_2 = total_sum_2 + recentered_sum_2
                    nUsed_2_orbit = nUsed_2_orbit + nTerms_2_module                                                      
                                                            
            # NORMALIZATION
            total_sum_1 = total_sum_1 / nUsed_1_orbit  
            total_sum_2 = total_sum_2 / nUsed_2_orbit  
            
            # BG SUBTRACTION
            bg_1 = imageSums_utilities.calculateBackground_noImg(total_sum_1)
            bgSubtracted_total_sum_1 = total_sum_1 - bg_1
            
            bg_2 = imageSums_utilities.calculateBackground_noImg(total_sum_2)
            bgSubtracted_total_sum_2 = total_sum_2 - bg_2

            ### INTEGRATE ###
            q = 2 * numpy.pi / orbit.resolution
            sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
            sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
            
            integratedIntensity_ellipse_1 = imageSums_utilities.integrate_ellipse(bgSubtracted_total_sum_1, sigma_x, sigma_y, 
                                                                                  ellipse_multiplicative_factor, expansion_factor=precision_factor)
            integratedIntensity_ellipse_1 = integratedIntensity_ellipse_1/((10**precision_factor)**2)
            
            integratedIntensity_ellipse_2 = imageSums_utilities.integrate_ellipse(bgSubtracted_total_sum_2, sigma_x, sigma_y, 
                                                                                  ellipse_multiplicative_factor, expansion_factor=precision_factor)
            integratedIntensity_ellipse_2 = integratedIntensity_ellipse_2/((10**precision_factor)**2)
    
            fOpen.write('%4d%4d%15.4f%15.4f%4d%4d\n'%(h_label, k_label, integratedIntensity_ellipse_1, integratedIntensity_ellipse_2, nUsed_1_orbit, nUsed_2_orbit))
                                                                      
        fOpen.close()
    

if __name__ == "__main__":
    print "\n**** CALLING imageSums_CChalf ****"
    imageSums_CChalf(sys.argv[1:])