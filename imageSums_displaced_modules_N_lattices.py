# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle
import getopt
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import makeOrbits
import imageSums_utilities
import buildReciprocalLattice
import calculatePredictedPattern

           
def imageSums(myArguments):
    str1 = '--selectedRun <selectedRun> --resolutionLimit <resolutionLimit>'
    str2 = '--halfWidth <halfWidth> --imagesDirectoryName <imagesDirectoryName> --geometryFile <geometryFile>'
    str3 = '--nominalCell <nominalCell> --hmax <hmax> --kmax <kmax>'
    str4 = '--tiltAngle <tiltAngle> --wavelength <wavelength> --detectorDistance  <detectorDistance> --pixelSize <pixelSize>'
    str5 = '--precision_factor <precision_factor> --N_lattices <N_lattices>'
    
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
                                                                 "pixelSize=",
                                                                 "precision_factor=",
                                                                 "N_lattices="])
                                                                 
    except getopt.GetoptError:
        print 'python imageSums_displaced_modules_N_lattices.py %s %s %s %s %s'%(str1, str2, str3, str4, str5)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':            
            print 'python imageSums_displaced_modules_N_lattices.py %s %s %s %s %s'%(str1, str2, str3, str4, str5)
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
        elif option == "--precision_factor":
            precision_factor = int(value)
        elif option == "--N_lattices":
            N_lattices = int(value)
            
    print 'N_lattices = ', N_lattices
    
    truncated_halfWidth = halfWidth - 10
    
    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_%d_lattices'%(selectedRun, N_lattices)   
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
    
    # EXTRACT MODULE DISPLACEMENTS
    module_displacements_file = open('./Output_r%s/ModuleDisplacements/module_displacements.pkl'%selectedRun, 'rb')
    modules = pickle.load(module_displacements_file) # bottom, top, left, right, N, <Dx0>, sig_x0, <Dy0>, sig_y0
    module_displacements_file.close()
    
    # ADD MODULE LIMITS (USED LATER, TO DISCARD MODULES OUT OF ORBIT RESOLUTION)
    modules_with_d_range = []
    for module in modules:
        bottom_bound = int(module[0])
        top_bound = int(module[1])
        left_bound = int(module[2])
        right_bound = int(module[3])
        print '\n*** MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound
        
        # CALCULATE MODULE RESOLUTION RANGE 
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
        
        minDistance = min([d_left_bottom, d_left_top, d_right_bottom, d_right_top])  - 0.01  # m
        maxDistance = max([d_left_bottom, d_left_top, d_right_bottom, d_right_top])  + 0.01  # m
        print 'MODULE DISTANCE FROM DIRECT BEAM RANGE: [%f, %f]'%(minDistance, maxDistance)
        
        module_with_d_range = [module[0], module[1], module[2], module[3], module[4], module[5], module[6], module[7], module[8], minDistance, maxDistance]
        modules_with_d_range.append(module_with_d_range)
    
    # CALCULATE STANDARD PATTERN        
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(nominalCell, hmax, kmax, resolutionLimit)    
    wavevector = 2 * numpy.pi/wavelength
    trialInPlaneRotation = 0
    standardPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, wavevector, tiltAngle, detectorDistance, pixelSize)

    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit) # 220 orbits to 4 A
          
    # FOR EACH ORBIT, FOR EACH MODULE, DO SECTOR SUMS
    orbits_withSum = {}
    for orbit in orbits:
        
        # FIGURES FOLDER
        if orbit.resolution > 7.0:
            orbit_figureFolder = '%s/low_res'%outputFolder                                      
        if 6.0 < orbit.resolution <= 7.0:
            orbit_figureFolder = '%s/high_res_7_6'%outputFolder            
        if 5.0 < orbit.resolution <= 6.0:
            orbit_figureFolder = '%s/high_res_6_5'%outputFolder            
        if 4.0 < orbit.resolution <= 5.0:
            orbit_figureFolder = '%s/high_res_5_4'%outputFolder
                  
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
        
        # EXTRACT ORBIT RADIUS ON DETECTOR
        spotDistance = numpy.nan
        for predictedSpot in standardPattern:
            if predictedSpot[0] == h_label and predictedSpot[1] == k_label:
                spotDistance = float(predictedSpot[10]) * pixelSize # m
                               
        print 'ORBIT: ', orbit.orbitIndices, 'LABEL: ', orbit.label, 'RADIUS ON DETECTOR: ', spotDistance
        
        # ONE ORBIT
        nTot_orbit = 0
        nUsed_orbit = 0
        total_sum = numpy.zeros((2*truncated_halfWidth*(10**precision_factor), 2*truncated_halfWidth*(10**precision_factor))) # 300x300
        
        # LOOP ON ALL MODULES
        for module in modules_with_d_range:
            
            min_distance = float(module[9])
            max_distance = float(module[10])
            
            bottom_bound = int(module[0])
            top_bound = int(module[1])
            left_bound = int(module[2])
            right_bound = int(module[3])
            x0 = halfWidth + float(module[5])
            y0 = halfWidth + float(module[7])

            print 'MODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound 
            print '<Dx0> = %.2f, <Dy0> = %.2f'%(float(module[5]), float(module[7]))
            
            if spotDistance < min_distance or spotDistance > max_distance:
                print 'UNUSED MODULE [%.2f, %.2f] (ORBIT RADIUS: %.2f m)'%(min_distance, max_distance, spotDistance)
                continue
            
            module_sum = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 50x50
            module_sum = numpy.matrix(module_sum)

            nTerms_module = 0
            
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
                        
                        # SPOT COORDINATES IN UNASSEMBLED MATRIX
                        i = int(spot[5])
                        j = int(spot[6])
                        
                        # CHECK MODULE
                        if not (bottom_bound <= i < top_bound):
                            continue
                        if not (left_bound <= j < right_bound):
                            continue
                        
                        ### EXTRACT 50x50 SECTOR ###
                        left_edge   = j - halfWidth
                        right_edge  = j + halfWidth
                        bottom_edge = i - halfWidth
                        top_edge    = i + halfWidth
                        if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                            continue
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
                        module_sum = module_sum + spotMatrix_rotated # SUM ON SINGLE MODULE
                        
                        nTerms_module = nTerms_module + 1
                        nTot_orbit = nTot_orbit + 1
            
            print 'nTerms_module: ', nTerms_module
            if nTerms_module == 0:
                continue
                                                                                                                
            recentered_sum = imageSums_utilities.recenter(module_sum, x0, y0, precision_factor, truncated_halfWidth) # 50x50 -> 300x300, NO BG SUB, NO NORMALIZATION
            total_sum = total_sum + recentered_sum
            nUsed_orbit = nUsed_orbit + nTerms_module     

        # CHECK N TERMS
        if nUsed_orbit != nTot_orbit:
            print 'PROBLEM: nUsed_orbit DIFFERENT FROM nTOt_orbit.'                                                   
                                                        
        # NORMALIZATION
        total_sum = total_sum / nUsed_orbit  
        
        # BG SUBTRACTION
        bg, n_bgPixels = imageSums_utilities.calculateBackground_nBg(total_sum)
        bgSubtracted_total_sum = total_sum - bg
        
        # GENERATE UPDATED ORBITS DICTIONARY
        orbit.nTerms = nUsed_orbit
        orbit.n_bgPixels = n_bgPixels
        orbit.bgSubtracted_total_sum = bgSubtracted_total_sum
        orbits_withSum['%d_%d'%(h_label, k_label)] = orbit
        
        ### SINGLE PLOT ###
        myFigureObject = matplotlib.pyplot.figure()
        myAxesImageObject = matplotlib.pyplot.imshow(bgSubtracted_total_sum, origin='lower', interpolation='nearest')
        matplotlib.pyplot.title('Orbit: %d %d'%(h_label, k_label))
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)        
        matplotlib.pyplot.savefig('%s/modules_sum_%d_%d_nTot_%d.png'
                                  %(orbit_figureFolder, h_label, k_label, nTot_orbit), dpi = 2*96 )
        matplotlib.pyplot.close()  
        
        ### SUM INTEGRATION ###
        integratedIntensity = imageSums_utilities.integrate(bgSubtracted_total_sum, expansion_factor=precision_factor)     
        integratedIntensity = integratedIntensity/((10**precision_factor)**2)
        
        # GAUSS FIT
        refined_sigma_x, refined_sigma_y, \
        refined_x0, refined_y0, \
        refined_amplitude, gauss_integral, \
        data, data_fitted = imageSums_utilities.do_gaussFit(bgSubtracted_total_sum)
        
        if not numpy.isnan(gauss_integral):
            refined_sigma_x = refined_sigma_x/(10**precision_factor)
            refined_sigma_y = refined_sigma_y/(10**precision_factor)
            gauss_integral = gauss_integral/((10**precision_factor)**2)
        
        print 'INTENSITY: ', integratedIntensity, gauss_integral
        
        ### PLOT GAUSS FIT ###
        
        if not numpy.isnan(gauss_integral):
            myFigureObject = matplotlib.pyplot.figure()
            myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*truncated_halfWidth*(10**precision_factor), 2*truncated_halfWidth*(10**precision_factor)), 
                                                         origin='lower', interpolation='nearest')
            try:
                matplotlib.pyplot.gca().contour(numpy.linspace(0, 2*truncated_halfWidth*(10**precision_factor)-1, 2*truncated_halfWidth*(10**precision_factor)), 
                                                numpy.linspace(0, 2*truncated_halfWidth*(10**precision_factor)-1, 2*truncated_halfWidth*(10**precision_factor)), 
                                                data_fitted.reshape(2*truncated_halfWidth*(10**precision_factor), 2*truncated_halfWidth*(10**precision_factor)), 4, colors='w')
            except:
                print 'PROBLEM DRAWING CONTOURS'
            
            matplotlib.pyplot.title('Orbit: %d %d Gauss integral: %.1f counts\nSig_x %.2f Sig_y %.2f'
                                     %(h_label, k_label, gauss_integral, refined_sigma_x, refined_sigma_y))
            myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
            matplotlib.pyplot.savefig('%s/gauss_module_sum_%d_%d_nTot_%d.png'
                                      %(orbit_figureFolder, h_label, k_label, nTot_orbit), dpi = 2*96)
            matplotlib.pyplot.close() 
        
        fOpen.write('%4d%4d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n'%(h_label, k_label, 
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