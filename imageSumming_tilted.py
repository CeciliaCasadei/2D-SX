# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle
import sys
import getopt
import matplotlib

import makeOrbits
import imageSums_utilities
import LP


            
def sumTilted_Function(myArguments):
    
    tiltAngles = {}
    tiltAngles['5']  = ['0198', '0199']
    tiltAngles['15'] = ['0195', '0196', '0197']
    tiltAngles['20'] = ['0200', '0201']

    # DEFAULTS:
    binStep = 0.01
    binSize_factor = 3
        
    str1 = '--binStep <binStep> --binSize_factor <binSize_factor>'
    str2 = '--resolutionLimit <resolutionLimit> --halfWidth <halfWidth>'
    str3 = '--imagesDirectoryName <imagesDirectoryName> --geometryFile <geometryFile>'
    str4 = '--nominalCell <nominalCell> --wavelength <wavelength>'
    str5 = '--tiltAngle <tiltAngle> --detectorDistance <detectorDistance>'
        
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["binStep=",
                                                                 "binSize_factor=", 
                                                                 "resolutionLimit=", 
                                                                 "halfWidth=", 
                                                                 "imagesDirectoryName=",
                                                                 "geometryFile=",
                                                                 "nominalCell=",
                                                                 "wavelength=", # 1.48 A
                                                                 "tiltAngle=",
                                                                 "detectorDistance=" # 0.285 m
                                                                 ])
    except getopt.GetoptError:
        print ('Error Usage: python imageSumming_tilted.py %s %s %s %s %s'
               %(str1, str2, str3, str4, str5))
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print ('Usage: python imageSumming_tilted.py %s %s %s %s %s'
                   %(str1, str2, str3, str4, str5))
            sys.exit()
        elif option == "--binStep":
            binStep = float(value)
        elif option == "--binSize_factor":
            binSize_factor = float(value)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--imagesDirectoryName":
            imagesDirectoryName = value
        elif option == "--geometryFile":
            geometryFile = value   
        elif option == "--nominalCell":
            cellSize = float(value)
        elif option == "--wavelength":
            wavelength = float(value)
        elif option == "--tiltAngle":
            tiltAngle_deg = value
        elif option == "--detectorDistance":
            detectorDistance = float(value)
            
    runs = tiltAngles['%s'%tiltAngle_deg]
    tiltAngle = float(tiltAngle_deg) * numpy.pi / 180
            
    waveVector = 2 * numpy.pi / wavelength     
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],
                                          [0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I  
    
    q_binWidth = (float(binSize_factor)/2)*binStep       
    
    # FOLDERS  
    outputFolder = './Output_imageSum_tilt_%s'%tiltAngle_deg  
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
    
    # EXTRACT GEOMETRY
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    
    # EXTRACT MODULE DISPLACEMENTS (bottom, top, left, right, N, <Dx0>, sig_Dx0, <Dy0>, sig_Dy0)
    module_displacements_file = open('./Output_r0127/ModuleDisplacements/module_displacements.pkl', 
                                     'rb')
    modules = pickle.load(module_displacements_file) 
    module_displacements_file.close()

    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
        #if orbit_label[0] == 2 and orbit_label[1] == 5:
            rodIndices.append(orbit)         
    print '%d orbits'%len(rodIndices)
            
    # LOOP ON ORBITS
    dictionaryList = []
    for orbit in rodIndices:

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
               
        reciprocalVector = [h_label, k_label]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        
        sinAlpha = 1
        qMax = (waveVector*(
                numpy.cos(tiltAngle)-numpy.sqrt(
                numpy.cos(tiltAngle)**2 - 
                (q_2D/waveVector)**2 - 2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha
                )
                ))
        print h_label, k_label, qMax        
        
        qLimits = []
        for n in range(0, 10000):
            qLimit = n * binStep
            if qLimit <= qMax:
                qLimits.append(qLimit)
            else:
                break     
        for i in range(1, len(qLimits)):
            qLimits.append(-qLimits[i])            
        qLimits.sort()    
        
        for q_item in range(0, len(qLimits)-1):
            middleQ = qLimits[q_item]
            leftQ  = middleQ - q_binWidth
            rightQ = middleQ + q_binWidth
           
            # SUM (FIXED h, k, qRod)
            sumMatrix_dictionary = {}
            sumMatrix = numpy.zeros(shape = (2*halfWidth, 2*halfWidth)) # 30x30 or 50x50
            sumMatrix = numpy.matrix(sumMatrix)    
            nTerms_sum = 0
           
            for run in runs:
                print run
                # PREPARE LATTICES TO IMAGES MATCHING
                lattices = joblib.load('./Output_runMergingVsModel/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                           %(run, run)) # 629 lattices  
                lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'
                                      %(run, run, run), 'r')
                lattices_names = list(lattices_names)    
                images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'
                                    %(run, run))
                images_names = list(images_names)
                # LOOP ON ALL PROCESSED LATTICES (FIXED ORBIT, FIXED MODULE)
                for index in range(0, len(lattices_names)):
                    lattice_name = lattices_names[index]
                    imageNumber = int(lattice_name[80:84])
                    latticeMatrix = lattices[index]                    # h_t k_t q_rod I_scaled flag i_unassembled j_unassembled scale
                    image_name = images_names[imageNumber-1].split()[1]
                    
                    if latticeMatrix[0, 4] == 0:                       # Check transformation and scaling flag
                        continue
                   
                    # LOAD UNASSEMBLED IMAGE
                    unassembledDataFile = h5py.File('%s/r%s-images/data1/%s'
                                                     %(imagesDirectoryName, run, image_name), 'r')
                    unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
                    unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####  (1480, 1552)
                    nRows = unassembledData.shape[0]
                    nColumns = unassembledData.shape[1]
                    
                    # IN THE CURRENT LATTICE, LOOP ON SPOTS BELONGING TO THE CURRENT ORBIT  
                    for spot in latticeMatrix:
                        h = int(spot[0])
                        k = int(spot[1])
                        qRod = float(spot[2])
                        if ((h == h_1 and k == k_1 and leftQ <= qRod <= rightQ) or 
                            (h == h_2 and k == k_2 and leftQ <= qRod <= rightQ) or 
                            (h == h_3 and k == k_3 and leftQ <= qRod <= rightQ) or
                            (-h == h_1 and -k == k_1 and leftQ <= -qRod <= rightQ) or 
                            (-h == h_2 and -k == k_2 and leftQ <= -qRod <= rightQ) or 
                            (-h == h_3 and -k == k_3 and leftQ <= -qRod <= rightQ)):
                                
                            i = int(spot[5])
                            j = int(spot[6])
                            
                            left_edge   = j - halfWidth
                            right_edge  = j + halfWidth
                            bottom_edge = i - halfWidth
                            top_edge    = i + halfWidth
                            if (left_edge < 0 or right_edge > nColumns or 
                                bottom_edge < 0 or top_edge > nRows):
                                continue
                            
                            ### EXTRACT SECTOR ###
                            spotMatrix = unassembledData[bottom_edge:top_edge, 
                                                         left_edge:right_edge]  # 30x30 or 50x50
                            
                            ### SCALE ###
                            spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                            spotMatrix = spot[7] * spotMatrix
                            
                            ### LP FACTOR ###
                            xDetector = xGeometry_np[i, j]
                            yDetector = yGeometry_np[i, j]
                            LPfactor = LP.get_LP(xDetector, 
                                                 yDetector,
                                                 q_2D,
                                                 qRod,
                                                 waveVector,
                                                 detectorDistance,
                                                 tiltAngle)        
                            
                            print LPfactor
                            spotMatrix = LPfactor * spotMatrix
                            
                            ### DETERMINE AZIMUTH ON DETECTOR ###
                            detectorAzimuth = imageSums_utilities.calculate_detectorAzimuth(xGeometry_np, 
                                                                                            yGeometry_np, 
                                                                                            i, 
                                                                                            j)
                                    
                            ### DETERMINE MODULE ROTATION ANGLE ###
                            moduleRotation = imageSums_utilities.calculate_moduleRotation(xGeometry_np, 
                                                                                          yGeometry_np, 
                                                                                          i, 
                                                                                          j)
                             
                            ### ROTATION ANGLE ###
                            rotationAngle = - detectorAzimuth + moduleRotation
                            rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                            
                            ### ROTATE ###  ### CLOCKWISE !!! ###
                            spotMatrix_rotated = imageSums_utilities.clockWiseRotation(spotMatrix, 
                                                                                       rotationAngle)
                            
                            ### GEOMETRY CORRECTION ###
                            module_Dx0 = numpy.nan
                            module_Dy0 = numpy.nan
                            for module_line in modules:
                                bottomBound = module_line[0]
                                topBound = module_line[1]
                                leftBound = module_line[2]
                                rightBound = module_line[3]
                                if (bottomBound <= i <= topBound and 
                                    leftBound <= j <= rightBound):
                                    module_Dx0 = module_line[5]
                                    module_Dy0 = module_line[7]
                
                            if (numpy.isnan(module_Dx0) or numpy.isnan(module_Dy0)):
                                print 'PROBLEM'
                    
                            recentered_spotMatrix = imageSums_utilities.translate(spotMatrix_rotated, 
                                                                                         module_Dx0, 
                                                                                         module_Dy0) # NO BG SUB, NO NORMALIZATION, DETECTOR COUNTS
                                        
                            ### SUM ###
                            sumMatrix = sumMatrix + recentered_spotMatrix # SUM ON SINGLE ORBIT, Qrod BIN                               
                            nTerms_sum = nTerms_sum + 1
                                            
            print 'h k q nTerms: ', h_label, k_label, middleQ, nTerms_sum
            if not nTerms_sum == 0:    
                bg, n_bgPixels = imageSums_utilities.calculateBackground_nBg(sumMatrix)
                bgSubtracted_total_sum = sumMatrix - bg
                bgSubtracted_total_sum = bgSubtracted_total_sum / nTerms_sum
        
                 
                myFigureObject = matplotlib.pyplot.figure()
                myAxesImageObject = matplotlib.pyplot.imshow(bgSubtracted_total_sum, 
                                                             origin='lower', 
                                                             interpolation='nearest')
                matplotlib.pyplot.title('Orbit: h = %d k = %d qRod = %.3f'
                                         %(h_label, k_label, middleQ))
                myFigureObject.colorbar(myAxesImageObject, 
                                        pad=0.01, 
                                        fraction=0.0471, 
                                        shrink=1.00, 
                                        aspect=20)        
                matplotlib.pyplot.savefig('%s/modules_sum_%d_%d_Q_%d_nTot_%d.png'
                                           %(outputFolder, h_label, k_label, q_item, nTerms_sum), 
                                           dpi = 2*96 )
                matplotlib.pyplot.close()  
                                           
                ####
                #### SAVE RESULTS
                #### 
                sumMatrix_dictionary['h'] = h_label
                sumMatrix_dictionary['k'] = k_label
                sumMatrix_dictionary['qRod'] = middleQ
                sumMatrix_dictionary['nTerms'] = nTerms_sum
                sumMatrix_dictionary['sumMatrix'] = bgSubtracted_total_sum
                dictionaryList.append(sumMatrix_dictionary)
                                                        
    sums_file = open('%s/sumMatrix_dictionary_list_tilt_%s.pkl'%(outputFolder, tiltAngle_deg), 'wb')
    pickle.dump(dictionaryList, sums_file)
    sums_file.close()

if __name__ == "__main__":
    print "\n**** CALLING imageSumming_tilted ****"
    sumTilted_Function(sys.argv[1:])   