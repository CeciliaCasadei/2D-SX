# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import h5py

import runSetup_settings
import discontinuityMask
import simulate_resolution
import LP
import buildMasks_module 



def calculate_nBg(boxMatrix, maskBoxMatrix):   

    # Returns n of pixels originally used to model the background                                
    
    gridStep = 6
    subBoxWidth = 6
    lowFluctuationThreshold = 2 # 2.5 in MATLAB
    
    xGrid = range(6, boxMatrix.shape[1], gridStep)   # 6 12 18 ... 90 (max 15 values)  
    yGrid = range(6, boxMatrix.shape[0], gridStep)   # 6 12 18 ... 90 (max 15 values)  
    
    stdDevMatrix = numpy.zeros((len(yGrid), len(xGrid)))   # 15x15 or smaller
    avgsMatrix   = numpy.zeros((len(yGrid), len(xGrid)))   # 15x15 or smaller
    stdDevVector = []
    
    boxMatrixNrows    = boxMatrix.shape[0]
    boxMatrixNcolumns = boxMatrix.shape[1]
    
    xGridItemIdx = 0        
    for xGridItem in xGrid:
        yGridItemIdx = 0            
        for yGridItem in yGrid:
            stdDevMatrix[yGridItemIdx, xGridItemIdx] = 1000
            avgsMatrix[yGridItemIdx, xGridItemIdx] = -1
            
            xLeft_subBox  = max([xGridItem - subBoxWidth, 0])
            xRight_subBox = min([xGridItem + subBoxWidth, boxMatrixNcolumns])
            yDown_subBox  = max(yGridItem - subBoxWidth, 0)
            yUp_subBox    = min(yGridItem + subBoxWidth, boxMatrixNrows)
            
            subBoxMatrix     = boxMatrix[yDown_subBox:yUp_subBox, 
                                         xLeft_subBox:xRight_subBox]           # 12x12 or smaller
            subMaskBoxMatrix = maskBoxMatrix[yDown_subBox:yUp_subBox, 
                                             xLeft_subBox:xRight_subBox]       # 12x12 or smaller
            myValidIndices = numpy.argwhere(subMaskBoxMatrix == 1)             # Exclude neighbouring modules
            nValid = myValidIndices.shape[0]
            
            if (nValid > 1 and 
                nValid > (subMaskBoxMatrix.shape[0] * subMaskBoxMatrix.shape[1] / 2)):                                        
                avg = numpy.sum(subBoxMatrix) / nValid              
                mySquare = numpy.square(subBoxMatrix)
                mySum = numpy.sum(mySquare)                                    # print mySum.dtype ---> float64
                myVariance = (mySum / nValid) - (avg ** 2)
                if myVariance > 0:
                    stdDev = numpy.sqrt(myVariance)
                    stdDevVector.append(stdDev)
                    stdDevMatrix[yGridItemIdx, xGridItemIdx] = stdDev
                    avgsMatrix[yGridItemIdx, xGridItemIdx]   = avg                    
                else:
                    print "BG SUBTRACTION PROBLEM: VARIANCE %.18f"%myVariance                    
            yGridItemIdx = yGridItemIdx + 1
        xGridItemIdx = xGridItemIdx + 1
    localNoise = numpy.percentile(stdDevVector, 15)
    
    lowFluctuationIndices \
    = numpy.argwhere(stdDevMatrix <= lowFluctuationThreshold * localNoise)
    
    sectorSize   = boxMatrix.shape[0] * boxMatrix.shape[1]
    n_gridPoints = len(xGrid) * len(yGrid)
    n_fitPoints  = len(lowFluctuationIndices)
    if n_fitPoints > n_gridPoints:
        print 'PROBLEM!!!'
    n_bgPixels = int(sectorSize * ( float(n_fitPoints) / n_gridPoints))  
                                                          
    return n_bgPixels


    

def calculate_individualErrorsFunction(myArguments):
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    boxWidth = 48
    nCountsPerPhoton = 26
    
    string_1 = '--imagesDirectoryPath <imagesDirectoryPath>'
    string_2 = '--geometryFile <geometryFile>'
    string_3 = '--nominalCell <nominalCell>'
    string_4 = '--wavelength <wavelength>'
    string_5 = '--detectorDistance <detectorDistance>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["imagesDirectoryPath=",
                                               "geometryFile=",
                                               "nominalCell=",
                                               "wavelength=",
                                               "detectorDistance="])
    except getopt.GetoptError:        
        print 'Usage: python calculate_individualErrors.py %s %s %s %s %s'%(string_1,
                                                                            string_2,
                                                                            string_3,
                                                                            string_4,
                                                                            string_5)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_individualErrors.py %s %s %s %s %s'%(string_1,
                                                                                string_2,
                                                                                string_3,
                                                                                string_4,
                                                                                string_5)
            sys.exit()
        elif option == "--imagesDirectoryPath":
            imagesDirectoryPath = value
        elif option == "--geometryFile":
            geometryFile = value
        elif option == "--nominalCell":
            cellSize = float(value)
        elif option == "--wavelength":
            wavelength = float(value)
        elif option == "--detectorDistance":
            detectorDistance = float(value)
    
    baseFolder = './Output_runMergingVsModel'
    waveVector = 2 * numpy.pi / wavelength
    tiltAngle_dict = runSetup_settings.tiltAngles
    
    ### EXTRACT GEOMETRY ###
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
      
    for runNumber in runNumbers:
        
        # SAVE RESULTS
        lattices_list_error = []
        
        # GET TILT ANGLE (Needed in LP factor calculation)
        tiltAngle_deg = tiltAngle_dict['%s'%runNumber]
        tiltAngle = float(tiltAngle_deg) * numpy.pi / 180
        
        # GET NORMALIZED SCALES (VS MODEL)
        folder = '%s/transformAndScaleToModel_r%s'%(baseFolder, runNumber)
        normalizedScales = joblib.load('%s/r%s-scales/r%s-normalizedScales.jbl'
                                        %(folder, runNumber, runNumber))
        
        # GET LATTICES (AFTER TRANSFORMATION VS MODEL)
        folder_t = '%s/spotsMatricesList-Transformed-r%s'%(folder, runNumber)
        latticesList = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                                    %(folder_t, runNumber))
        
        # GET LATTICES (AFTER SCALING VS MODEL)
        folder_s = '%s/spotsMatricesList-Scaled-r%s'%(baseFolder, runNumber)
        latticesList_s = joblib.load('%s/r%s_scaledSpotsMatricesList.jbl'
                                      %(folder_s, runNumber))
        
        # GET IMAGE PATH
        lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'
                                      %(runNumber, runNumber, runNumber), 'r')
        lattices_names = list(lattices_names)    
        images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'
                                    %(runNumber, runNumber))
        images_names = list(images_names)
        
        print '****************'  
        print 'RUN ', runNumber
        print len(normalizedScales), ' LATTICES' 
        print len(latticesList), ' LATTICES'  
        print len(lattices_names), ' LATTICES'  
        
        for lattice in range(0, len(normalizedScales)):
            print lattice
            
            # SAVE RESULTS
            lattice_matrix_error = []
            
            # EXTRACT NORMALIZED SCALE
            normalizedScale = float(normalizedScales[lattice])
            
            # EXTRACT LATTICE MATRIX
            latticeMatrix = latticesList[lattice]  # h (transformed vs model)
                                                   # k (transformed vs model)
                                                   # qRod
                                                   # I (Photons, LP-corrected, un-Scaled)
                                                   # flag
                                                   # i (unassembled matrix)
                                                   # j (unassembled matrix)
            
            # EXTRACT SCALED LATTICE MATRIX
            latticeMatrix_s = latticesList_s[lattice]  # h (transformed vs model)
                                                       # k (transformed vs model)
                                                       # qRod
                                                       # I (Photons, LP-corrected, Scaled)
                                                       # flag
                                                       # i (unassembled matrix)
                                                       # j (unassembled matrix)
                                                 
            # EXTRACT RAW IMAGE
            lattice_name = lattices_names[lattice]
            imageNumber = int(lattice_name[80:84])                    
            image_name = images_names[imageNumber-1].split()[1]
                   
            # LOAD UNASSEMBLED IMAGE
            unassembledDataFile = h5py.File('%s/r%s-images/data1/%s'
                                             %(imagesDirectoryPath, runNumber, image_name), 'r')
            unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
            unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####  (1480, 1552)
            unassembledDataNrows    = unassembledData.shape[0]
            unassembledDataNcolumns = unassembledData.shape[1]
            
            for n in range(latticeMatrix.shape[0]):
                spot   = latticeMatrix[n]
                spot_s = latticeMatrix_s[n]
                                
                # EXTRACT INTENSITY
                I    = spot[3]    # (Photons, LP-corrected, un-Scaled)
                I_S  = spot_s[3]  # (Photons, LP-corrected, Scaled)
                
                flag = spot_s[4]  # Flag is zero when transform vs model or
                                  #                   scale vs model
                                  # did not succeed.
                                  # Flag can be one but with nan intensity
                                  # when integration area was not on a single detector module.
                                  # (integrateSpot.py)
                                  
                if flag == 1:
                    if (abs(I_S - normalizedScale*I) > 0.001) :
                        raise Exception('PROBLEM')
                        
                # SPOT COOS IN UNASSEMBLED MATRIX
                i = int(spot[5])
                j = int(spot[6])
                
                # EXTRACT SECTOR
                xLeft  = max([j - boxWidth, 0])
                xRight = min([j + boxWidth, unassembledDataNcolumns])
                yDown  = max([i - boxWidth, 0])
                yUp    = min([i + boxWidth, unassembledDataNrows])   
                
                boxMatrix = unassembledData[yDown:yUp, xLeft:xRight]   
                            
                ### BUILD DETECTOR MODULE MASK ###
                maskBoxMatrix = numpy.ones((2*boxWidth, 2*boxWidth), 
                                            dtype=numpy.int)                             # 96x96 
                maskBorder_left, \
                maskBorder_down, \
                maskBoxMatrix = discontinuityMask.discontinuityMask(xGeometry_np, 
                                                                    yGeometry_np, 
                                                                    maskBoxMatrix, 
                                                                    i, 
                                                                    j, 
                                                                    boxWidth)  # maskBoxMatrix 96x96 or smaller
                maskedBoxMatrix = numpy.multiply(boxMatrix, maskBoxMatrix)     # 96x96 or smaller, zero in neighbouring modules
                                                                               # ORIGINAL DETECTOR COUNTS
            
                # GET N OF PIXELS ORIGINALLY USED TO CALCULATE BACKGROUND
                n_bgPxls = calculate_nBg(maskedBoxMatrix, maskBoxMatrix)   
                
                # GET N OF PIXELS ORIGINALLY USED TO INTEGRATE
                integrationMask, \
                ringMask = buildMasks_module.buildMasks(maskedBoxMatrix, 
                                                        1, 
                                                        5, 
                                                        5)        
                n_integrationPxls = integrationMask.sum()  
                
                n_ring = ringMask.sum()   
            
                # CONVERT SINGLE IMG IN N OF PH - NO SCALE - NO LP
                spotMatrix_ph = boxMatrix / nCountsPerPhoton
                
                # LP FACTOR
                h    = spot[0]
                k    = spot[1]
                qRod = spot[2]
                if h != spot_s[0] or k != spot_s[1] or qRod != spot_s[2]:
                    raise Exception('ERROR')
                d_2D = simulate_resolution.resolution(cellSize, h, k, 0)
                q_2D = 2*numpy.pi/d_2D
                
                xDetector = xGeometry_np[i, j]
                yDetector = yGeometry_np[i, j]
                
                LP_factor = LP.get_LP(xDetector, 
                                      yDetector,
                                      q_2D,
                                      qRod,
                                      waveVector,
                                      detectorDistance,
                                      tiltAngle)                
                                       
                ### CALCULATE variance IN BG REGION (PHOTONS, LP-CORRECTED, NO SCALE)
                spotMatrix_ph = LP_factor * spotMatrix_ph
                integratedRing = numpy.multiply(ringMask, spotMatrix_ph).sum()
                avg_integratedRing = integratedRing / n_ring
                ring_intensities = numpy.multiply(ringMask, spotMatrix_ph)                    
                ring_intensities = ring_intensities.flatten().T                    
            
                sum_sq = 0
                m = 0
                for ring_pxl in ring_intensities:
                    
                    if ring_pxl != 0: # In ring region, however pixels with zero value in ring region are excluded                     
                        m = m+1
                        diff = ring_pxl - avg_integratedRing
                        sq = diff**2
                        sum_sq = sum_sq + sq
                
                # Consider pixels with zero value in ring region
                sum_sq = sum_sq + (n_ring - m)*(avg_integratedRing**2)
                
                var_ring = sum_sq/n_ring  
                
                
                # CALCULATE ERROR ON UNSCALED I
                r = n_integrationPxls/n_bgPxls
                dI_squared = I + (r + 1)*n_integrationPxls*var_ring
                if dI_squared < 0: # I can be negative ....
                    dI_squared = (r + 1)*n_integrationPxls*var_ring
                dI = numpy.sqrt(dI_squared)
                           
                # CALCULATE ERROR ON SCALED I
                dI_S = normalizedScale * dI                
                #print flag, I, I_S, dI_S
                
                lattice_matrix_error.append([spot_s[0],                     # h
                                             spot_s[1],                     # k
                                             spot_s[2],                     # qRod
                                             spot_s[3],                     # I, LP_corrected, Photons, scaled
                                             spot_s[4],                     # flag (0 if T vs model or S vs model failed)
                                             spot_s[5],                     # i
                                             spot_s[6],                     # j
                                             dI_S])                         # Error on scaled I
            lattice_matrix_error = numpy.asarray(lattice_matrix_error)      # 2D matrix with 8 columns
            lattices_list_error.append(lattice_matrix_error)
        
        print "************"
        print len(lattices_list_error)
        print "************"
        joblib.dump(lattices_list_error, '%s/r%s_scaledSpotsMatricesList_error.jbl'
                                          %(folder_s, runNumber))
                

     
if __name__ == "__main__":
    print "\n**** CALLING calculate_individualErrors ****"
    calculate_individualErrorsFunction(sys.argv[1:])    