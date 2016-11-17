# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle
import scipy.interpolate

import makeOrbits
import imageSums_utilities

    

def resolutionCutoff_singleImage():
    
    # PARAMETERS
    selectedRun = '0127'
    resolutionLimit = 4.0
    halfWidth = 15
    multiplicative_factor = 2.5
    nCountsPerPhoton = 26
    
    # FOLDERS
    outputFolder = './Output_resolutionCutoff_singleImage'    
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
        
    # EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
    fRead = open('./Output_imageSums_sigmaFits/sigmaXCurveParameters.pkl', 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()
    
    fRead = open('./Output_imageSums_sigmaFits/sigmaYCurveParameters.pkl', 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close()
               
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
    
    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                   # 220 orbits to 4 A
    
    # LISTS STORING RESULTS
    ratios = []
    qs = []
    
    # LOOP ON ORBITS
    for orbit in orbits:
        
        print 'ORBIT: ', orbit.orbitIndices, 'LABEL: ', orbit.label
        
        indices_1 = orbit.orbitIndices[0]
        indices_2 = orbit.orbitIndices[1]
        indices_3 = orbit.orbitIndices[2]
        h_1 = indices_1[0]
        k_1 = indices_1[1]
        h_2 = indices_2[0]
        k_2 = indices_2[1]
        h_3 = indices_3[0]
        k_3 = indices_3[1]
        
        # LOOP ON ALL PROCESSED LATTICES 
        for index in range(0, len(lattices_names)):
            lattice_name = lattices_names[index]
            runNumber = lattice_name[10:14]
            imageNumber = int(lattice_name[80:84])
            latticeNumber = lattice_name[92]
            latticeMatrix = lattices[index]                    # h_t k_t q_rod I_scaled flag i_unassembled j_unassembled scale
            image_name = images_names[imageNumber-1].split()[1]
            
            if latticeMatrix[0, 4] == 0:                       # Check transformation and scaling flag
                continue
            
            print runNumber, imageNumber, latticeNumber, image_name
            
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
                    left_edge   = j - halfWidth
                    right_edge  = j + halfWidth
                    bottom_edge = i - halfWidth
                    top_edge    = i + halfWidth
                    if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                        continue
                    
                    ### EXTRACT 30x30 SECTOR ###
                    spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]
                    
                    ### SCALE ###
                    spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                    spotMatrix = spot[7] * spotMatrix
                    
                    ### DETERMINE AZIMUTH ON DETECTOR ###
                    xDetector = xGeometry_np[i, j]
                    yDetector = yGeometry_np[i, j]
                    if xDetector != 0:
                        detectorAzimuth = numpy.arctan(yDetector/xDetector)
                        if xDetector < 0 and yDetector > 0:
                            detectorAzimuth = detectorAzimuth + numpy.pi
                        if xDetector < 0 and yDetector < 0:
                            detectorAzimuth = detectorAzimuth - numpy.pi
                    else:
                        if yDetector > 0:
                            detectorAzimuth = + numpy.pi /2
                        else:
                            detectorAzimuth = - numpy.pi /2
                            
                    ### DETERMINE MODULE ROTATION ANGLE ###
                    deltaXgeo = float(xGeometry_np[i, j+1] - xGeometry_np[i, j])
                    deltaYgeo = float(yGeometry_np[i, j+1] - yGeometry_np[i, j])
    
                    if deltaXgeo != 0:
                        moduleRotation = numpy.arctan(deltaYgeo/deltaXgeo) # From module to lab frame
                        if deltaXgeo < 0 and deltaYgeo > 0:
                            moduleRotation = moduleRotation + numpy.pi
                        elif deltaXgeo < 0 and deltaYgeo < 0:
                            moduleRotation = moduleRotation - numpy.pi
                    else:
                        if deltaYgeo > 0:
                            moduleRotation = numpy.pi / 2
                        else:
                            moduleRotation = -numpy.pi / 2
                     
                    ### ROTATION ANGLE ###
                    rotationAngle = - detectorAzimuth + moduleRotation
                    rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                    ### ROTATE ###  ### CLOCKWISE !!! ###
                    x_windows_pix = range(-spotMatrix.shape[1]/2, +spotMatrix.shape[1]/2)  # -18, -17, ..., +17
                    y_windows_pix = range(-spotMatrix.shape[0]/2, +spotMatrix.shape[0]/2)  # -18, -17, ..., +17
    
                    [X_windows_pix, Y_windows_pix] = numpy.meshgrid(x_windows_pix, y_windows_pix)
                    X_windows_pix = numpy.asarray(X_windows_pix, dtype=numpy.float32)
                    Y_windows_pix = numpy.asarray(Y_windows_pix, dtype=numpy.float32)
                    
                    X_windows_pix_rotated = numpy.cos(rotationAngle)*X_windows_pix - numpy.sin(rotationAngle)*Y_windows_pix
                    Y_windows_pix_rotated = numpy.sin(rotationAngle)*X_windows_pix + numpy.cos(rotationAngle)*Y_windows_pix
                    
                    f = scipy.interpolate.interp2d(x_windows_pix, y_windows_pix, spotMatrix, kind='linear')
                    spotMatrix_rotated = numpy.zeros(spotMatrix.shape)
                    
                    for columnIndex in range(0, spotMatrix_rotated.shape[1]):
                        for rowIndex in range(0, spotMatrix_rotated.shape[0]):
                            rotated_x = X_windows_pix_rotated[rowIndex, columnIndex]
                            rotated_y = Y_windows_pix_rotated[rowIndex, columnIndex]
                            rotated_f = f(rotated_x, rotated_y)
                            spotMatrix_rotated[rowIndex, columnIndex] = rotated_f                       
                                       
                    ### CALCULATE BACKGROUND ###
                    background_rotated = imageSums_utilities.calculateBackground_noImg(spotMatrix_rotated)
                       
                    ### BG SUBTRACTION ###
                    bgSubtracted_rotated = spotMatrix_rotated - background_rotated
                    
                    # I/SIG(I) CALCULATION
                    q = 2 * numpy.pi / orbit.resolution
                    sigma_x = imageSums_utilities.quadratic(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
                    sigma_y = imageSums_utilities.line_plus_sigmoid(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2], sigmaYCurveParameters[3], sigmaYCurveParameters[4])    
                               
                    # N COUNTS TO N PHOTONS CONVERSION
                    bgSubtracted_rotated =  bgSubtracted_rotated/nCountsPerPhoton
                    
                    # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
                    integrationMask = numpy.zeros((bgSubtracted_rotated.shape))   
                    ringMask = numpy.zeros((bgSubtracted_rotated.shape)) 
                    
                    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
                    x_axis = multiplicative_factor * sigma_x
                    y_axis = multiplicative_factor * sigma_y
                    centerX = integrationMask.shape[1] / 2
                    centerY = integrationMask.shape[0] / 2
                    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
                    
                    integrationMask[numpy.where(distance < 1)] = 1
                    ringMask[numpy.where(distance > 1.5)] = 1
                    ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
                    
                    n_peak = integrationMask.sum()
                    n_bg = ringMask.sum()
                    
                    integratedIntensity = numpy.multiply(integrationMask, bgSubtracted_rotated).sum()
                    avg_I = integratedIntensity/n_peak
                       
                    integratedBg = numpy.multiply(ringMask, bgSubtracted_rotated).sum()
                    avg_bg = integratedBg/n_bg
                    ring_bg = numpy.multiply(ringMask, bgSubtracted_rotated)
                    ring_bg = ring_bg.flatten().T
                    
                    sum_sq = 0
                    n = 0
                    for bg_pxl in ring_bg:
                        if bg_pxl != 0:
                            n = n+1
                            diff = bg_pxl - avg_bg
                            sq = diff**2
                            sum_sq = sum_sq + sq
                    
                    bg_fluctuations = numpy.sqrt(sum_sq/n_bg)
                    ratio = avg_I/bg_fluctuations
                
                    qs.append(q)
                    ratios.append(ratio)
    
    # SAVE LISTS
    qs_file = open('%s/qs.pkl'%outputFolder, 'wb')
    pickle.dump(qs, qs_file)
    qs_file.close()    
    
    ratios_file = open('%s/ratios.pkl'%outputFolder, 'wb')
    pickle.dump(ratios, ratios_file)
    ratios_file.close()    
    

if __name__ == "__main__":
    print "\n**** CALLING resolutionCutoff_singleImage ****"
    resolutionCutoff_singleImage()