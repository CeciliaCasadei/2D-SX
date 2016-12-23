# -*- coding: utf-8 -*-



"""
@author: casadei_c
xGeometry[i,j] is the x coordinate on detector (in m, wrt center) of intensity stored in unassembledData[i,j]  
yGeometry[i,j] is the y coordinate on detector (in m, wrt center) of intensity stored in unassembledData[i,j]  
"""



### STANDARD ###
import sys
import getopt
import pickle
import numpy
import math
import h5py
import os
import time
import joblib
import scipy.optimize
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
from matplotlib import rcParams



### PYTHON ###
import peakDetection
import buildReciprocalLattice
import backGroundSubtraction   # Cython plane subtraction in planeSubtraction.pyx
import calculatePredictedPattern



### CYTHON ###
import discontinuityMask
import unassembledMatching
import calculateLatticeError
import recalculateDistance



### CLASSES ###
import spotClass



def toBeMinimized(x, *p):
    unitCell = x[0]/4
    orientation = x[1]/285
    imageCenter = [x[2]/10, x[3]/10]
    
    highResLimit, waveVector, tiltAngle, detectorDistance, pixelSize, distanceThreshold, detectedPeaks = p
    
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(unitCell, 100, 100, highResLimit)
    predictedPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, orientation, 
                                                                                   waveVector, tiltAngle, 
                                                                                   detectorDistance, pixelSize)    
    predictedPattern = numpy.asarray(predictedPattern, dtype=numpy.float32)
    imageCenter = numpy.asarray(imageCenter, dtype=numpy.float32)   
    nMatchedPeaks, latticeError = calculateLatticeError.calculateMatrixElement(imageCenter, predictedPattern, 
                                                         detectedPeaks, pixelSize, 
                                                         distanceThreshold)
    return latticeError



def processing(myArguments):
    
    ### DEFAULTS ###
    runNumber = ''
    bgSubtractionMethod = 'plane'
    minimizationMethod = '4Dbf'
    fractionDetectedThreshold = 0
    lowResLimit = 50.0
    highResLimit = 7.1
    nCountsPerPhoton = 26
    integrationRadius = 5
    geometryFile= ''
    imageFolder = ''
    imageSelection = ''
    latticeSelection = ''
    
    # READ INPUTS
    string1 = 'Usage: python processingAndIntegration.py --runNumber <runNumber> --bgSubtractionMethod <bgSubtractionMethod>'
    string2 = ' --minimizationMethod <minimizationMethod> --fractionDetectedThreshold <fractionDetectedThreshold> --lowResLimit <lowResLimit> --highResLimit <highResLimit>'
    string3 = ' --nCountsPerPhoton <nCountsPerPhoton> --integrationRadius <integrationRadius> --geometryFile <geometryFile>'
    string4 = ' --imageFolder <imageFolder> --imageSelection <imageSelection> --latticeSelection <latticeSelection>'    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", 
                                                                 "bgSubtractionMethod=", "minimizationMethod=", 
                                                                 "fractionDetectedThreshold=",
                                                                 "lowResLimit=", "highResLimit=", 
                                                                 "nCountsPerPhoton=", "integrationRadius=",
                                                                 "geometryFile=", "imageFolder=",
                                                                 "imageSelection=", "latticeSelection="])
    except getopt.GetoptError:
        print string1 + string2 + string3 + string4
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':            
            print string1 + string2 + string3 + string4
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--bgSubtractionMethod":
            bgSubtractionMethod = value
        elif option == "--minimizationMethod":
            minimizationMethod = value
        elif option == "--fractionDetectedThreshold":
            fractionDetectedThreshold = float(value)            
        elif option == "--lowResLimit":
            lowResLimit = float(value)    
        elif option == "--highResLimit":
            highResLimit = float(value)   
        elif option == "--nCountsPerPhoton":
            nCountsPerPhoton = int(value)  
        elif option == "--integrationRadius":
            integrationRadius = int(value)
        elif option == "--geometryFile":
            geometryFile = value
        elif option == "--imageFolder":
            imageFolder = value
        elif option == "--imageSelection":
            imageSelection = value
        elif option == "--latticeSelection":
            latticeSelection = int(value)
                  
    ### OTHER PARAMETERS ###
    nIterations = 3
    nSizeRefSteps = 15
    nOrientationRefSteps = 15
    nCenterRefSteps = 15
    widthSizeRefSteps = [0.0012, 0.0008, 0.0004]                               #[0.0016, 0.0010, 0.0006]
    widthOrientationRefSteps = [0.12, 0.08, 0.04]                              #[0.16, 0.10, 0.06]
    widthOrientationRefSteps = numpy.asarray(widthOrientationRefSteps)
    widthOrientationRefSteps = widthOrientationRefSteps/180*numpy.pi
    widthCenterRefSteps = [0.06, 0.04, 0.02]                                   #[0.10, 0.06, 0.04]
    distanceThreshold = 8
    distanceThresholds = [5, 4, 3]
    boxWidth = 48
     
    ### FOLDERS ###
    processingFolder = './Output_r%s/UnassembledImageProcessing'%runNumber        
    if not os.path.exists(processingFolder):
        os.mkdir(processingFolder)
    processingFiguresFolder = '%s/BgSubtractionAndPeakSearchPlots'%processingFolder
    
    ### EXTRACT GEOMETRY ###
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    
    ### EXTRACT LATTICES ###    
    latticeObjectsPath = './Output_r%s/OrientationAndCellRefinement/r%s_refinedLatticesDictionary.pkl'%(runNumber, runNumber)
    fLattices = open(latticeObjectsPath, 'rb')
    latticesDictionary = pickle.load(fLattices)
    fLattices.close()
    
    
    ### STORE REFINED LATTICE SIZE ###
    refinedLatticeSizes = []
    
    ### LOOP ON LATTICES ###
    for myKey, myLattice in latticesDictionary.items():
        processFlag = 0
        # IF ONE IMAGE / LATTICE WAS SELECTED, SWITCH FIGURE PLOTS ON.
        if not imageSelection == '' and not latticeSelection == '':
            
            if not os.path.exists(processingFiguresFolder):
                os.mkdir(processingFiguresFolder)
            if myLattice.runNumber == runNumber and myLattice.imageNumber == imageSelection and myLattice.latticeNumberInImage == latticeSelection:
                processFlag = 1
                bgPlanePlotFlag = 1
                singleSpotFigureFlag = 1
                integratedPeaksFigureFlag = 1
        else:
            if myLattice.runNumber == runNumber:
                processFlag = 1
                bgPlanePlotFlag = 0
                singleSpotFigureFlag = 0
                integratedPeaksFigureFlag = 0
                 
        if processFlag == 1:
            print "PROCESSING - Run %s Image %s Lattice %s"%(myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage)
            
            ### SET IMAGE CENTER ###
            imageCenter = [0, 0]
            imageCenter = numpy.asarray(imageCenter, dtype=numpy.float32)    
            myLattice.imageCenter = imageCenter
            
            ### EXTRACT LATTICE ATTRIBUTES ###
            tiltAngle = float(myLattice.tiltAngle)/180*numpy.pi
            waveVector = 2*numpy.pi/myLattice.wavelength
            detectorDistance = myLattice.detectorDistance
            pixelSize = myLattice.pixelSize
            
            ### INITIALIZE LISTS TO FOLLOW LATTICE REFINEMENT BEHAVIOUR ###
            refinedCellSizes = []
            refinedCellSizes.append(myLattice.refinedCellSize)
            refinedLatticeOrientations = []
            refinedLatticeOrientations.append(myLattice.refinedInPlaneOrientation)
            refinedCenterXs = []
            refinedCenterXs.append(myLattice.imageCenter[0])
            refinedCenterYs = []
            refinedCenterYs.append(myLattice.imageCenter[1])
            nDetectedAndMatchedPeaks = []
            nDetectedAndMatchedPeaks.append(myLattice.nMatchedPeaks)
            avgLatticeErrors = []
            avgLatticeErrors.append(myLattice.avgLatticeError)
            
            ### EXTRACT MEASUREMENT RAW/UNASSEMBLED DATA ####
            unassembledDataFile = h5py.File('%s/%s'%(imageFolder, myLattice.fileName), 'r')
            unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
            unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####
            unassembledDataNrows = unassembledData.shape[0]
            unassembledDataNcolumns = unassembledData.shape[1]

            ### LOGGING ###
            fOpen = open('%s/Processing_r%s_img%s_lattice%s_%s.txt'%(processingFolder,
                                                                     myLattice.runNumber, 
                                                                     myLattice.imageNumber, 
                                                                     myLattice.latticeNumberInImage,
                                                                     minimizationMethod), 'w') 
            fOpen.write('Run: %s \nImage: %s \nLattice: %s\n'%(myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage))
            fOpen.write('START: \nLattice orientation: %.5f \nCell size: %.5f \n'%(myLattice.refinedInPlaneOrientation, myLattice.refinedCellSize)) 
            fOpen.write('Center x: %.6f\n'%myLattice.imageCenter[0])
            fOpen.write('Center y: %.6f\n'%myLattice.imageCenter[1])
            fOpen.write('\n***** Connected peaks detection *****\n\n')  
            fOpen.write('    n    h    k    Predicted x    Predicted y    Local noise     Detected x     Detected y     Detected I    Distance (pxls)\n\n')
            
            ### LOOP ON PREDICTED PEAKS ###
            print 'Raw sectors extraction, background subtraction and connected peaks detection.\n'
            startTime = time.time()
            
            nPredictedPeak = 0  # Spot ID number     
            
            detectedPeaks  = [] # LIST ONLY SPOTS INITIALLY FOUND CLOSE TO PREDICTION (DISTANCE BELOW distanceThreshold) - TO BE USED IN 4D REFINEMENT
            spotDictionary = {} # COLLECT ***ALL*** SPOTS IN RESOLUTION RANGE
            
            refinedPredictedPattern = myLattice.refinedPredictedPattern        # h k qx qy dmin q 
                                                                               # azimuth rotated_azimuth detector_azimuth 
                                                                               # diffraction_angle detector_radius qrod LPfactor       
            for predictedPeak in refinedPredictedPattern:
                if predictedPeak[4] <= lowResLimit and predictedPeak[4] >= highResLimit:
                    nPredictedPeak = nPredictedPeak + 1
                    
                    h_peak = predictedPeak[0]
                    k_peak = predictedPeak[1]
                    diffractionSpot = spotClass.diffractionSpot(nPredictedPeak, h_peak, k_peak)
                    
                    ### PREDICTED SPOT POSITION (INITIAL CENTER IN 0,0) ###
                    Xpredicted = pixelSize * predictedPeak[10] * numpy.cos(predictedPeak[8])  ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                    Ypredicted = pixelSize * predictedPeak[10] * numpy.sin(predictedPeak[8])  ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                    
                    ### MATCH PREDICTED SPOT POSITION TO POSITION IN UNASSEMBLED MATRIX VIA GEOMETRY ###
                    xIndices = numpy.argwhere(abs(xGeometry - Xpredicted) <= 0.000055)
                    xIndices = numpy.asarray(xIndices, dtype = numpy.int16)
                    yIndices = numpy.argwhere(abs(yGeometry - Ypredicted) <= 0.000055)
                    yIndices = numpy.asarray(yIndices, dtype = numpy.int16)
                    successFlag, i_good, j_good = unassembledMatching.unassembledMatching(xIndices, yIndices)
                    # xGeometry[i_good,j_good] ~ Xpredicted AND yGeometry[i_good,j_good] ~ Ypredicted -> (Xpredicted, Ypredicted) correspond to unassembledData[i_good, j_good]  
                    
                    ### EXTRACT UNASSEMBLED DATA SECTOR ###
                    if successFlag == 0:
                        fOpen.write('%5d%5d%5d \tPredicted position was not found in geometry.\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1]))  
                        spotDictionary['Spot%d'%nPredictedPeak] = diffractionSpot 
                    else:
                        xLeft  = max([j_good - boxWidth, 0])
                        xRight = min([j_good + boxWidth, unassembledDataNcolumns])
                        yDown  = max([i_good - boxWidth, 0])
                        yUp    = min([i_good + boxWidth, unassembledDataNrows])            
                        boxMatrix = unassembledData[yDown:yUp, xLeft:xRight]   # 96x96 or smaller unassembled data sector
                        
                        ### BUILD DETECTOR MODULE MASK ###
                        maskBoxMatrix = numpy.ones((2*boxWidth, 2*boxWidth), dtype=numpy.int)                             # 96x96 
                        maskBorder_left, maskBorder_down, maskBoxMatrix = discontinuityMask.discontinuityMask(xGeometry_np, yGeometry_np, 
                                                                                                              maskBoxMatrix, 
                                                                                                              i_good, j_good, 
                                                                                                              boxWidth)   # maskBoxMatrix 96x96 or smaller
                        maskedBoxMatrix = numpy.multiply(boxMatrix, maskBoxMatrix)                                        # 96x96 or smaller, zero in neighbouring modules
                        
                        ### BACKGROUND PLANE SUBTRACTION ###
                        bgSubtractionDictionary = backGroundSubtraction.backGroundSubtractionFunction(maskedBoxMatrix, maskBoxMatrix, bgSubtractionMethod, 
                                                                                                      myLattice.fileName, myLattice.runNumber, myLattice.imageNumber, 
                                                                                                      myLattice.latticeNumberInImage, nPredictedPeak, 
                                                                                                      bgPlanePlotFlag, processingFiguresFolder)            
                        localNoise            = bgSubtractionDictionary['localNoise']
                        bgSubtractedBoxMatrix = bgSubtractionDictionary['bgSubtractedBoxMatrix']      # 96x96 or smaller
                        bgSubtractedBoxMatrix = numpy.multiply(bgSubtractedBoxMatrix, maskBoxMatrix)  # 96x96 or smaller, zero in neighbouring modules
                        
                        ### DETERMINATION OF MODULE ROTATION ANGLE ###
                        if j_good+1 < xGeometry.shape[1]:                
                            deltaXgeo = float(xGeometry[i_good, j_good+1] - xGeometry[i_good, j_good])
                            deltaYgeo = float(yGeometry[i_good, j_good+1] - yGeometry[i_good, j_good])
                        else:
                            deltaXgeo = float(xGeometry[i_good, j_good] - xGeometry[i_good, j_good-1])
                            deltaYgeo = float(yGeometry[i_good, j_good] - yGeometry[i_good, j_good-1])
                        if deltaXgeo != 0:
                            moduleRotation = numpy.arctan(deltaYgeo/deltaXgeo) # Rotation required to convert module coos to real space
                            if deltaXgeo < 0 and deltaYgeo > 0:
                                moduleRotation = moduleRotation + numpy.pi
                            elif deltaXgeo < 0 and deltaYgeo < 0:
                                moduleRotation = moduleRotation - numpy.pi
                        else:
                            if deltaYgeo > 0:
                                moduleRotation = numpy.pi / 2
                            else:
                                moduleRotation = -numpy.pi / 2
                            
                        ### DETERMINATION OF ORIGIN IN BOX (BELONGING TO CURRENT MODULE) ###
                        labX_Origin = xGeometry[yDown+maskBorder_down, xLeft+maskBorder_left]
                        labY_Origin = yGeometry[yDown+maskBorder_down, xLeft+maskBorder_left] 
                        
                        ### CONNECTED PEAKS DETECTION ###
                        peakDetectionDictionary = peakDetection.peakDetectionFunction(bgSubtractedBoxMatrix, 
                                                                                      localNoise, pixelSize,
                                                                                      moduleRotation, 
                                                                                      labX_Origin, labY_Origin,
                                                                                      maskBorder_left, maskBorder_down,
                                                                                      Xpredicted, Ypredicted)                    
                        peakDetectionMask             = peakDetectionDictionary['peakDetectionMask']
                        xCenterOfMass_vector          = peakDetectionDictionary['xCenterOfMass_vector'] # In lab frame, in m and relative to detector center, to be compared to xPredicted
                        yCenterOfMass_vector          = peakDetectionDictionary['yCenterOfMass_vector'] # In lab frame, in m and relative to detector center, to be compared to yPredicted
                        integratedIntensity_vector    = peakDetectionDictionary['integratedIntensity_vector']
                        distanceFromPrediction_vector = peakDetectionDictionary['distanceFromPrediction_vector']
                        
                        diffractionSpot.setBoxMatrix(boxMatrix)
                        diffractionSpot.setMaskBoxMatrix(maskBoxMatrix)
                        diffractionSpot.setBgSubtractedBoxMatrix(bgSubtractedBoxMatrix)
                        diffractionSpot.setPeakDetectionMask(peakDetectionMask)
                        diffractionSpot.setBoxLimits(xLeft, yDown)                 
                        diffractionSpot.setBoxIndices(i_good-yDown, j_good-xLeft)
                        
                        nConnectedPeaks = len(integratedIntensity_vector)
                        if nConnectedPeaks == 0:
                            spotDictionary['Spot%d'%nPredictedPeak] = diffractionSpot 
                            fOpen.write('%5d%5d%5d \tNo detected peaks.\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1]))
                        else:
                            minIdx = distanceFromPrediction_vector.index(min(distanceFromPrediction_vector))
                            if distanceFromPrediction_vector[minIdx] > distanceThreshold:
                                spotDictionary['Spot%d'%nPredictedPeak] = diffractionSpot 
                                fOpen.write('%5d%5d%5d \tDetected peak(s) too far from prediction.\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1]))
                            else:                                                            ### KEEP THIS EXPERIMENTAL PEAK (relaxed threshold) ###
                                detectedPeak = []
                                xCoM = xCenterOfMass_vector[minIdx]
                                yCoM = yCenterOfMass_vector[minIdx]
                                detectedPeak.append(h_peak)                                  # 0 ---> h
                                detectedPeak.append(k_peak)                                  # 1 ---> k
                                detectedPeak.append(xCoM)                                    # 2 ---> x CoM in m, in lab frame (to be compared with Xpredicted)
                                detectedPeak.append(yCoM)                                    # 3 ---> y CoM in m, in lab frame (to be compared with Ypredicted)
                                detectedPeak.append(integratedIntensity_vector[minIdx])      # 4 ---> Integrated intensity
                                detectedPeak.append(distanceFromPrediction_vector[minIdx])   # 5 ---> Discrepancy (pxls) between prediction and observation (center of mass)
                                detectedPeaks.append(detectedPeak)
                                                         
                                diffractionSpot.setCoM(xCoM, yCoM)
                                diffractionSpot.setConnectedI(integratedIntensity_vector[minIdx])
                                diffractionSpot.setDistanceFromPrediction(distanceFromPrediction_vector[minIdx])
                                
                                spotDictionary['Spot%d'%nPredictedPeak] = diffractionSpot 
                                fOpen.write('%5d%5d%5d%15.8f%15.8f%15.2f%15.8f%15.8f%15.2f%19.2f\n\n'%(nPredictedPeak, h_peak, k_peak, 
                                                                                                       Xpredicted, Ypredicted, localNoise, xCoM, yCoM,                                                                                                                    
                                                                                                       integratedIntensity_vector[minIdx], 
                                                                                                       distanceFromPrediction_vector[minIdx])) 
                                                     
            detectedPeaks = numpy.asarray(detectedPeaks, dtype=numpy.float32)
            nDictionaryItems = 0
            for key, value in spotDictionary.items():
                nDictionaryItems = nDictionaryItems + 1
            
            if nDictionaryItems != nPredictedPeak:
                print 'PROBLEM!'
            nDetected = detectedPeaks.shape[0]   
            fOpen.write('Detected peaks (initial distance from prediction below %d pxls): %d/%d'%(distanceThreshold, nDetected, nDictionaryItems))
            
            endTime = time.time() - startTime
            print 'It took %.2f s\n'%endTime
            
            nDetectedThreshold = nDictionaryItems * fractionDetectedThreshold
            print 'N predicted: %d'%nDictionaryItems
            print 'Threshold: %d'%nDetectedThreshold
            print 'N detected: %d'%nDetected
            if nDetected <= nDetectedThreshold:
                print 'N of detected peaks is below threshold.'
                continue
            else:
                print 'N of detected peaks is above threshold.'
                
            ### 4D REFINEMENT ###
            print 'Refinement started. Method: %s\n'%minimizationMethod
            fOpen.write('\nREFINEMENT METHOD; %s \n'%minimizationMethod)
            refinementStart = time.time()
            
            calculateMatrixElement = calculateLatticeError.calculateMatrixElement
            nIteration = 0
            while nIteration < nIterations:
                print 'N iteration: %d'%nIteration 
                
                ### 4-DIMENSIONS BRUTE-FORCE ###
                if minimizationMethod == '4Dbf':
                    
                    ### SET SEARCH INTERVALS ###
                    nSizeRefStepsOverTwo = nSizeRefSteps / 2
                    sizeRefinementRescalings = numpy.linspace(1-nSizeRefStepsOverTwo*widthSizeRefSteps[nIteration], 
                                                              1+nSizeRefStepsOverTwo*widthSizeRefSteps[nIteration], 
                                                              nSizeRefSteps)
                    sizeRefinementSteps = myLattice.refinedCellSize * sizeRefinementRescalings
                
                    nOrientationRefStepsOverTwo = nOrientationRefSteps / 2
                    orientationRefinementSteps = numpy.linspace(-nOrientationRefStepsOverTwo*widthOrientationRefSteps[nIteration],
                                                                +nOrientationRefStepsOverTwo*widthOrientationRefSteps[nIteration],
                                                                nOrientationRefSteps)
                    orientationRefinementSteps = myLattice.refinedInPlaneOrientation + orientationRefinementSteps
                    
                    centerX_0 = myLattice.imageCenter[0]
                    centerY_0 = myLattice.imageCenter[1]
                    nCenterRefStepsOverTwo = nCenterRefSteps / 2
                    centerXs = centerX_0 + numpy.linspace(-nCenterRefStepsOverTwo*widthCenterRefSteps[nIteration],
                                                          +nCenterRefStepsOverTwo*widthCenterRefSteps[nIteration],
                                                          nCenterRefSteps)
                    centerYs = centerY_0 + numpy.linspace(-nCenterRefStepsOverTwo*widthCenterRefSteps[nIteration],
                                                          +nCenterRefStepsOverTwo*widthCenterRefSteps[nIteration],
                                                          nCenterRefSteps)
                    
                    ### BUILD 4D LATTICE-ERROR MATRIX ###
                    latticeErrorMatrix = numpy.zeros((nOrientationRefSteps, nSizeRefSteps, nCenterRefSteps, nCenterRefSteps))
                    
                    nCalls = 0
                    nSizeStep = 0
                    for sizeRefinementStep in sizeRefinementSteps: # (21) trial unit cell sizes
                        nSizeStep = nSizeStep + 1
                        reciprocalLatticeData = buildReciprocalLattice.buildReciprocalLatticeFunction(sizeRefinementStep, 100, 100, highResLimit)
                        
                        nRotationStep = 0       
                        for trialInPlaneRotation in orientationRefinementSteps: # (21) trial in-plane orientations
                            nRotationStep = nRotationStep + 1
                            
                            predictedPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLatticeData, trialInPlaneRotation, waveVector, tiltAngle, 
                                                                                                           myLattice.detectorDistance, myLattice.pixelSize)    
                            predictedPattern = numpy.asarray(predictedPattern, dtype=numpy.float32)
                    
                            centerXstep = 0    
                            for centerX in centerXs: 
                                centerXstep = centerXstep + 1
                                
                                centerYstep = 0
                                for centerY in centerYs: 
                                    centerYstep = centerYstep + 1
                                    
                                    trialImageCenter = [centerX, centerY]
                                    trialImageCenter = numpy.asarray(trialImageCenter, dtype=numpy.float32)           
                                    nMatchedPeaks, latticeError = calculateMatrixElement(trialImageCenter, predictedPattern, 
                                                                                         detectedPeaks, pixelSize, 
                                                                                         distanceThresholds[nIteration])
                                    nCalls = nCalls + 1
                                    latticeErrorMatrix[nRotationStep-1, nSizeStep-1, centerXstep-1, centerYstep-1] = latticeError
                     
                    print 'N calls %d'%nCalls               
                    i, j, k, l = numpy.unravel_index(latticeErrorMatrix.argmin(), latticeErrorMatrix.shape)
                    
                    ### RESULTS ###
                    refinedOrientation = orientationRefinementSteps[i]
                    refinedCellSize = sizeRefinementSteps[j]
                    refinedCenterX = centerXs[k]
                    refinedCenterY = centerYs[l]
                    
                else:
                    xStart = [myLattice.refinedCellSize*4, myLattice.refinedInPlaneOrientation*285, myLattice.imageCenter[0]*10, myLattice.imageCenter[1]*10]
                    fixedParas = (highResLimit, waveVector, tiltAngle, detectorDistance, pixelSize, distanceThresholds[nIteration], detectedPeaks)
                    result = scipy.optimize.minimize(toBeMinimized, xStart, args = fixedParas, method = minimizationMethod, 
                                                     options = {'xtol': 0.001, 'ftol': 0.001}) # Defaults: 'xtol': 0.0001, 'ftol': 0.0001
                    print result.message
                    
                    ### RESULTS ###
                    refinedOrientation = result.x[1]/285
                    refinedCellSize = result.x[0]/4
                    refinedCenterX = result.x[2]/10
                    refinedCenterY = result.x[3]/10
                
                ### CALCULATE & STORE ITERATION RESULTS ###
                reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(refinedCellSize, 100, 100, highResLimit)
                predictedPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, refinedOrientation, 
                                                                                               waveVector, tiltAngle, 
                                                                                               detectorDistance, pixelSize)    
                predictedPattern = numpy.asarray(predictedPattern, dtype=numpy.float32)
                imageCenter = [refinedCenterX, refinedCenterY]
                imageCenter = numpy.asarray(imageCenter, dtype=numpy.float32)
                nMatchedPeaks, latticeError = calculateMatrixElement(imageCenter, predictedPattern, 
                                                                     detectedPeaks, pixelSize, 
                                                                     distanceThresholds[nIteration])
                avgMinError = latticeError / nMatchedPeaks
                
                myLattice.setRefinedPattern(predictedPattern)
                myLattice.setRefinedInPlaneOrientation(refinedOrientation)
                myLattice.setRefinedCellSize(refinedCellSize)
                myLattice.imageCenter = imageCenter
                myLattice.setLatticeError(avgMinError)     
                            
                refinedLatticeOrientations.append(myLattice.refinedInPlaneOrientation)
                refinedCellSizes.append(myLattice.refinedCellSize)
                refinedCenterXs.append(myLattice.imageCenter[0])
                refinedCenterYs.append(myLattice.imageCenter[1])
                nDetectedAndMatchedPeaks.append(nMatchedPeaks)
                avgLatticeErrors.append(myLattice.avgLatticeError)
                
                ### LOG ITERATION RESULTS ###
                fOpen.write('\n*****')
                fOpen.write('\nIteration: %d'%nIteration)
                fOpen.write('\nLattice orientation: %.5f\n'%myLattice.refinedInPlaneOrientation)
                fOpen.write('Cell size: %.5f\n'%myLattice.refinedCellSize)
                fOpen.write('Center x: %.6f\n'%myLattice.imageCenter[0])
                fOpen.write('Center y: %.6f'%myLattice.imageCenter[1])
                
                ### AFTER EACH ITERATION, UPDATE DISTANCE FROM PREDICTION IN detectedPeaks TABLE ###
                detectedPeaks = recalculateDistance.recalculateDistanceFunction(imageCenter, predictedPattern, detectedPeaks, pixelSize)
                nIteration = nIteration + 1
                ### END 4D REFINEMENT ###
                                
            fOpen.close()
            refinementEnd = time.time() - refinementStart
            print 'Refinement took %.2f s\n'%refinementEnd  
            
            ### STORE REFINED LATTICE SIZE ###
            refinedLatticeSizes.append(myLattice.refinedCellSize)
            
            ### PLOT REFINEMENT BEHAVIOUR ###    
            myLattice.refinedLatticeOrientations = refinedLatticeOrientations
            myLattice.refinedCellSizes           = refinedCellSizes
            myLattice.refinedCenterXs            = refinedCenterXs
            myLattice.refinedCenterYs            = refinedCenterYs
            myLattice.nDetectedAndMatchedPeaks   = nDetectedAndMatchedPeaks
            myLattice.avgLatticeErrors           = avgLatticeErrors   
            
            myLattice.refinementBehaviourPlot(processingFolder, minimizationMethod)
            
            
            ### INTEGRATE AND PRODUCE SINGLE SPOT FIGURES ###
            integratedPeaks = []
            
            
            for mySpotKey, mySpot in spotDictionary.items():            
                h = mySpot.h
                k = mySpot.k
                
                successFlag = 0
                for predictedPeak in myLattice.refinedPredictedPattern:
                    if predictedPeak[0] == h and predictedPeak[1] == k:
                        refinedX = pixelSize * ( myLattice.imageCenter[0] + predictedPeak[10] * numpy.cos(predictedPeak[8]) )   ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                        refinedY = pixelSize * ( myLattice.imageCenter[1] + predictedPeak[10] * numpy.sin(predictedPeak[8]) )   ### POSITION IN m RELATIVE TO DIRECT BEAM ### 
                        xIndices = numpy.argwhere(abs(xGeometry - refinedX) <= 0.000055)
                        xIndices = numpy.asarray(xIndices, dtype = numpy.int16)
                        yIndices = numpy.argwhere(abs(yGeometry - refinedY) <= 0.000055)
                        yIndices = numpy.asarray(yIndices, dtype = numpy.int16)
                        successFlag, iFinal_global, jFinal_global = unassembledMatching.unassembledMatching(xIndices, yIndices)
                        qRod = predictedPeak[11]
                        LPfactor = predictedPeak[12]
                        break
                        
                if successFlag == 1:
                    if not hasattr(mySpot, 'yDown'):
                        ### ONE COULD LOOK FOR IT AGAIN ... ###
                        continue
                    else:                    
                        iFinal = iFinal_global-mySpot.yDown
                        jFinal = jFinal_global-mySpot.xLeft
                        if 0 <= iFinal < mySpot.peakDetectionMask.shape[0] and 0 <= jFinal < mySpot.peakDetectionMask.shape[0]:
                            mySpot.setFinalBoxIndices(iFinal, jFinal)
                            mySpot.nCountsPerPhoton = nCountsPerPhoton
                            mySpot.integrationRadius = integrationRadius
                            mySpot.integrateSpot()
                            
                            mySpot.qRod = qRod
                            correctedIntensity = LPfactor * mySpot.integratedIntensity
                            mySpot.correctedIntensity = correctedIntensity
                            
                            integratedPeak = []
                            integratedPeak.append(mySpot.n)
                            integratedPeak.append(mySpot.h)
                            integratedPeak.append(mySpot.k)
                            integratedPeak.append(mySpot.qRod)
                            integratedPeak.append(mySpot.integratedIntensity)
                            integratedPeak.append(mySpot.correctedIntensity)
                            
                            ####################################
                            integratedPeak.append(iFinal_global)
                            integratedPeak.append(jFinal_global)
                            ####################################
                            
                            integratedPeaks.append(integratedPeak)
                            spotDictionary['Spot%d'%mySpot.n] = mySpot  
                
                        if singleSpotFigureFlag == 1:
                            mySpot.spotProcessingPlot(myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage, boxWidth)
            
            ### SORTING ###
            integratedPeaks = numpy.asarray(integratedPeaks)
            n_data = integratedPeaks[:, 0]
            h_data = integratedPeaks[:, 1]
            k_data = integratedPeaks[:, 2]
            qRod_data = integratedPeaks[:, 3]
            I_data = integratedPeaks[:, 4]
            Ic_data = integratedPeaks[:, 5]
            
            ###################################
            iIndex_data = integratedPeaks[:, 6]
            jIndex_data = integratedPeaks[:, 7]
            ###################################
            
            for i in range(0, len(integratedPeaks)):
                for j in range(i+1, len(integratedPeaks)):
                    if n_data[i] > n_data[j]:
                        nTemp = n_data[i]
                        hTemp = h_data[i]
                        kTemp = k_data[i]
                        qRodTemp = qRod_data[i]
                        iTemp = I_data[i]
                        icTemp = Ic_data[i]
                        
                        ###########################
                        iIndexTemp = iIndex_data[i]
                        jIndexTemp = jIndex_data[i]
                        ###########################
                        
                        n_data[i] = n_data[j]
                        h_data[i] = h_data[j]
                        k_data[i] = k_data[j]
                        qRod_data[i] = qRod_data[j]
                        I_data[i] = I_data[j]
                        Ic_data[i] = Ic_data[j]
                        
                        ###############################
                        iIndex_data[i] = iIndex_data[j]
                        jIndex_data[i] = jIndex_data[j]
                        ###############################
                        
                        n_data[j] = nTemp
                        h_data[j] = hTemp
                        k_data[j] = kTemp
                        qRod_data[j] = qRodTemp
                        I_data[j] = iTemp  
                        Ic_data[j] = icTemp
                        
                        ###########################
                        iIndex_data[j] = iIndexTemp
                        jIndex_data[j] = jIndexTemp
                        ###########################
                        
            orderedIntegratedIntensities = numpy.zeros((len(integratedPeaks), 7))
            for i in range(0, len(integratedPeaks)):
                    orderedIntegratedIntensities[i,0]=h_data[i]
                    orderedIntegratedIntensities[i,1]=k_data[i]
                    orderedIntegratedIntensities[i,2]=qRod_data[i]
                    orderedIntegratedIntensities[i,3]=Ic_data[i]
                    orderedIntegratedIntensities[i,4]=1                        # flag
                    orderedIntegratedIntensities[i,5]=iIndex_data[i]
                    orderedIntegratedIntensities[i,6]=jIndex_data[i]
                    
            myLattice.orderedIntegratedIntensities = orderedIntegratedIntensities
            
            ### LOGGING ###
            integrationFile = open('%s/integration_r%s_img%s_lattice%s.txt'%(processingFolder, myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage), 'w')
            integrationFile.write('    n    h    k        qRod      I (ph)    I corrected\n\n')        
            for index in range(0, len(n_data)):
                integrationFile.write('%5d%5d%5d%12.5f%12.3f%15.3f\n'%(n_data[index], h_data[index], 
                                                                       k_data[index], qRod_data[index],
                                                                       I_data[index], Ic_data[index]))
            integrationFile.close()

            ### MINIMUM OUTPUT ###
            imageNumber = myLattice.imageNumber.zfill(4)
            joblib.dump(orderedIntegratedIntensities, '%s/OrderedIntegratedIntensities_r%s_Img%s_Lattice%s.jbl'
                                                       %(processingFolder, myLattice.runNumber, imageNumber, myLattice.latticeNumberInImage))    
                                                                                                            
            ### PAPER FIGURE (IMG N 74) ###
            ### NB Final orientation will be p -> LABEL INDICES NEED TO BE PERMUTED ###
            if integratedPeaksFigureFlag == 1:
                goodPredictions_x = []
                goodPredictions_y = []
                goodPredictions_labels = []
                badPredictions_x = []
                badPredictions_y = []
                badPredictions_labels = []
                zoomedSpots_x = []
                zoomedSpots_y = []
                zoomedSpots_labels = []
                                                           
                for spot in myLattice.refinedPredictedPattern:
                    h = spot[0]
                    k = spot[1]
                    label = r'%d, %d $\rightarrow$ %d, %d'%(h, k, k, h)
                    xPredicted = 874 + myLattice.refinedCenterXs[-1] + spot[10] * numpy.cos(spot[8]) 
                    yPredicted = 874 + myLattice.refinedCenterYs[-1] + spot[10] * numpy.sin(spot[8]) 
                    integrationFlag = 0
                    for mySpotKey, mySpot in spotDictionary.items():  
                        if mySpot.h == h and mySpot.k == k:
                            try:
                                I = mySpot.integratedIntensity
                                if not numpy.isnan(I):
                                    integrationFlag = 1
                            except:
                                integrationFlag = 0
                    if integrationFlag == 0:
                        badPredictions_x.append(xPredicted)
                        badPredictions_y.append(yPredicted)
                        badPredictions_labels.append(label)
                    else:
                        goodPredictions_x.append(xPredicted)
                        goodPredictions_y.append(yPredicted)
                        goodPredictions_labels.append(label)
                    if h == -13 and k == 11 or h == -11 and k == 13 or h == -13 and k == 2 or h == -11 and k == -2:
                        if h == -13 and k == 2:
                            label = r'-13, 2  $\equiv$ 2, 11 $\rightarrow$ 11, 2'
                        if h == -11 and k == 13:
                            label = r'-11, 13  $\equiv$ -2, -11 $\rightarrow$ -11, -2'
                        if h == -13 and k == 11:
                            label = r'-13, 11  $\equiv$ 11, 2 $\rightarrow$ 2, 11'
                        zoomedSpots_x.append(xPredicted)
                        zoomedSpots_y.append(yPredicted)
                        zoomedSpots_labels.append(label)
                        

                ### EXTRACT ASSEMBLED DATA FOR FINAL PLOTTING ###
                assembledData   = unassembledDataFile['/data/assembleddata0'] 
                assembledData   = numpy.asarray(assembledData, dtype=numpy.float32)         #### !!!!! ####
                
                ### INDEXED PREDICTED PATTERN FIGURE ###
                matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
                matplotlib.pyplot.title('%s'%myLattice.fileName, y=1.05)
                matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                
                matplotlib.pyplot.scatter(badPredictions_x, badPredictions_y, color='b', marker="o", linewidth='2', facecolors='none', s=600)
                for label, x, y in zip(badPredictions_labels, badPredictions_x, badPredictions_y):
                    matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                               textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                               bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
                matplotlib.pyplot.scatter(goodPredictions_x, goodPredictions_y,   color='r', marker="o", linewidth='2', facecolors='none', s=600)
                for label, x, y in zip(goodPredictions_labels, goodPredictions_x, goodPredictions_y):
                    matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                               textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                               bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
                
                matplotlib.pyplot.savefig("%s/r%s_Image_%s_Lattice_%s.png"%(processingFolder, runNumber, imageNumber, myLattice.latticeNumberInImage), fontsize = 20, dpi=96*4)
                matplotlib.pyplot.close()
                
                ### ZOOMED-IN FIGURE ###                
                matplotlib.pyplot.figure(figsize=(50,40), dpi=4*96, facecolor='w',frameon=True)
                rcParams['xtick.direction'] = 'in'
                rcParams['ytick.direction'] = 'in'
                matplotlib.pyplot.title('%s'%myLattice.fileName, y=1.05)
                ax1 = matplotlib.pyplot.subplot2grid((5,4), (0,0), rowspan=4, colspan=4)
                ax1.imshow(assembledData[250:1300, 700:1750], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                             
                badPredictions_x = numpy.asarray(badPredictions_x)
                badPredictions_x = badPredictions_x - 700
                goodPredictions_x = numpy.asarray(goodPredictions_x)
                goodPredictions_x = goodPredictions_x - 700     
                
                badPredictions_y = numpy.asarray(badPredictions_y)
                badPredictions_y = badPredictions_y - 250
                goodPredictions_y = numpy.asarray(goodPredictions_y)
                goodPredictions_y = goodPredictions_y - 250
                
                ax1.scatter(badPredictions_x,  badPredictions_y,  color='b', marker="o", linewidth='2', facecolors='none', s=600)
                ax1.scatter(goodPredictions_x, goodPredictions_y, color='r', marker="o", linewidth='2', facecolors='none', s=600)
                
                resolutionCircle = 7.0 #A
                myRadius = detectorDistance/pixelSize*math.tan(2*math.asin(myLattice.wavelength/(2*resolutionCircle)))
                circle = matplotlib.pyplot.Circle((874-700, 874-250),myRadius,linewidth=1.0, color='b',fill=False)
                ax1.add_artist(circle)
                
                ax1.set_xticklabels([])
                ax1.set_yticklabels([])
                ax1.set_xlim([0, 1050])
                ax1.set_ylim([0, 1050])
                
                ax2 = matplotlib.pyplot.subplot2grid((5,4), (4,0), colspan=1)
                ax2.imshow(assembledData[zoomedSpots_y[0]-21:zoomedSpots_y[0]+22, zoomedSpots_x[0]-21:zoomedSpots_x[0]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, color='b',fill=False)
                ax2.add_artist(circle)
                ax2.set_xticklabels([])
                ax2.set_yticklabels([])
                ax3 = matplotlib.pyplot.subplot2grid((5,4), (4,1), colspan=1)
                ax3.imshow(assembledData[zoomedSpots_y[1]-21:zoomedSpots_y[1]+22, zoomedSpots_x[1]-21:zoomedSpots_x[1]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, color='b',fill=False)
                ax3.add_artist(circle)
                ax3.set_xticklabels([])
                ax3.set_yticklabels([])
                ax4 = matplotlib.pyplot.subplot2grid((5,4), (4,2), colspan=1)
                ax4.imshow(assembledData[zoomedSpots_y[2]-21:zoomedSpots_y[2]+22, zoomedSpots_x[2]-21:zoomedSpots_x[2]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, color='b',fill=False)
                ax4.add_artist(circle)
                ax4.set_xticklabels([])
                ax4.set_yticklabels([])
                ax5 = matplotlib.pyplot.subplot2grid((5,4), (4,3), colspan=1)
                ax5.imshow(assembledData[zoomedSpots_y[3]-21:zoomedSpots_y[3]+22, zoomedSpots_x[3]-21:zoomedSpots_x[3]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, color='b',fill=False)
                ax5.add_artist(circle)
                ax5.set_xticklabels([])
                ax5.set_yticklabels([])
                           
                for i, ax in enumerate(matplotlib.pyplot.gcf().axes):
                    if i > 0:
                        ax.text(0.5, 0.5, "%s"%(zoomedSpots_labels[i-1]), horizontalalignment='left', verticalalignment='bottom', fontsize=28)          
                
                matplotlib.pyplot.savefig("%s/r%s_Image_%s_Lattice_%s_zoomed.png"%(processingFolder, runNumber, imageNumber, myLattice.latticeNumberInImage), fontsize = 20, dpi=96*4)
                matplotlib.pyplot.close()
                
    joblib.dump(refinedLatticeSizes, '%s/refinedLatticeSizes_r%s.jbl'%(processingFolder, runNumber))
                
if __name__ == "__main__":
    print "\n**** CALLING processing ****"
    processing(sys.argv[1:])