# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle
import sys
import getopt

import makeOrbits
import imageSums_utilities
import detectorModules

import buildReciprocalLattice
import calculatePredictedPattern
import partialSumClass


            
def calculate_partialSums_Function(myArguments):
    
    # DEFAULTS:
    selectedRun = '0127'
    resolutionLimit = 4.0
    halfWidth = 15
    
    str1 = '--selectedRun <selectedRun> --resolutionLimit <resolutionLimit>'
    str2 = '--halfWidth <halfWidth> --imagesDirectoryName <imagesDirectoryName> --geometryFile <geometryFile>'
    str3 = '--nominalCell <nominalCell> --hmax <hmax> --kmax <kmax>'
    str4 = '--tiltAngle <tiltAngle> --wavelength <wavelength> --detectorDistance  <detectorDistance> --pixelSize <pixelSize>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "resolutionLimit=", 
                                                                 "halfWidth=", 
                                                                 "imagesDirectoryName=",
                                                                 "geometryFile=",
                                                                 "nominalCell=",
                                                                 "hmax=", "kmax=",
                                                                 "tiltAngle=",
                                                                 "wavelength=",
                                                                 "detectorDistance=",
                                                                 "pixelSize="])
    except getopt.GetoptError:
        print 'Error Usage: python calculate_moduleDisplacements.py %s %s %s %s %s'%(str1, str2, str3, str4)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_moduleDisplacements.py %s %s %s %s %s'%(str1, str2, str3, str4)
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--imagesDirectoryName":
            imagesDirectoryName = value
        elif option == "--geometryFile":
            geometryFile = value
        elif option == "--nominalCell":
            nominalCell = float(value)
        elif option == "--hmax":
            hmax = int(value)
        elif option == "--kmax":
            kmax = int(value)
        elif option == "--tiltAngle":
            tiltAngle = float(value)
        elif option == "--wavelength":
            wavelength = float(value)
        elif option == "--detectorDistance":
            detectorDistance = float(value)
        elif option == "--pixelSize":
            pixelSize = float(value)
            
    # FOLDERS
    outputFolder = './Output_r%s/ModuleDisplacements'%selectedRun    
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
            
    # PREPARE LATTICES TO IMAGES MATCHING
    lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                           %(selectedRun, selectedRun, selectedRun)) # 629 lattices  
    lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(selectedRun, selectedRun, selectedRun), 'r')
    lattices_names = list(lattices_names)    
    images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(selectedRun, selectedRun))
    images_names = list(images_names)
    
    # EXTRACT GEOMETRY
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
    
    print 'xGeometry SHAPE: ', xGeometry_np.shape
    print 'yGeometry SHAPE: ', yGeometry_np.shape
    print 'MODULE ROW INTERFACES: ', row_interface
    print 'MODULE COLUMN INTERFACES: ', column_interface

    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                    # 220 orbits to 4 A
    print '%d orbits'%len(orbits)
        
    # CALCULATE STANDARD PATTERN    
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(nominalCell, hmax, kmax, resolutionLimit)
    tiltAngle = tiltAngle/180 * numpy.pi
    trialInPlaneRotation = 0
    wavevector = 2 * numpy.pi/wavelength
    standardPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, wavevector, tiltAngle, detectorDistance, pixelSize)
        
    # LOOP ON MODULES
    partialSums_list = []
    for module_row in range(0, len(row_interface)-1):
        for module_column in range(0, len(column_interface)-1):
            bottom_bound = row_interface[module_row]
            top_bound = row_interface[module_row+1]
            left_bound = column_interface[module_column]
            right_bound = column_interface[module_column+1]
                      
            print '\n*** MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound
            
            # CALCULATE MODULE RESOLUTION RANGE (speed up X4)
            border = 5 # Safe extraction from geometry.
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
            #print 'MODULE DISTANCE FROM DIRECT BEAM RANGE: [%f, %f]'%(minDistance, maxDistance)
            
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
                
                #print '\nOrbit %d, %d'%(h_label, k_label)       
                
                # CHECK ORBIT DETECTOR DISTANCE
                orbit_use_flag = 0
                
                spotDistance = numpy.nan
                for predictedSpot in standardPattern:
                    if predictedSpot[0] == h_label and predictedSpot[1] == k_label:
                        spotDistance = float(predictedSpot[10]) * pixelSize # m
                        if minDistance <= spotDistance <= maxDistance:
                            orbit_use_flag = 1
          
                if orbit_use_flag == 0:
                    print 'NOT using orbit %d, %d'%(h_label, k_label) 
                    #print 'Spot distance: %f out of [%f, %f]'%(spotDistance, minDistance, maxDistance)
                else:
                    #print 'Using orbit %d, %d'%(h_label, k_label)
                    #print 'Spot distance: %f in [%f, %f]'%(spotDistance, minDistance, maxDistance)
                    
                    # CALCULATE PARTIAL SUM (FIXED ORBIT, FIXED MODULE)
                    module_sum = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 30x30 or 50x50
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
                                
                                ### EXTRACT SECTOR ###
                                spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]  # 30x30 or 50x50
                                
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
                                                    
                    #print 'nTerms_module: ', nTerms_module_sum
                    if not nTerms_module_sum == 0:                                           
                        ####
                        #### SAVE OBJECT
                        #### ATTRIBUTES:
                        #### nTerms_module_sum
                        #### ORBIT LABEL
                        #### MODULE LABEL
                        #### module_sum (NOT BGSUB, NOT NORMALIZED, IN COUNTS)
                        ####
                        partialSumObject = partialSumClass.partialSum(h_label, k_label, 
                                                                      bottom_bound, top_bound, left_bound, right_bound, 
                                                                      nTerms_module_sum, module_sum)
                        partialSums_list.append(partialSumObject)
                                                        
    partialSums_file = open('%s/partialSums_list.pkl'%(outputFolder), 'wb')
    pickle.dump(partialSums_list, partialSums_file)
    partialSums_file.close()

if __name__ == "__main__":
    print "\n**** CALLING calculate_partialSums ****"
    calculate_partialSums_Function(sys.argv[1:])   