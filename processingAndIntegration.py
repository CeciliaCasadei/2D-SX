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
import h5py
import os
import time
import joblib
import scipy.optimize

### PYTHON ###
#import plotGeometry
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
    latticeSelection = 0
    
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
    widthSizeRefSteps = [0.0016, 0.0010, 0.0006]
    widthOrientationRefSteps = [0.16, 0.10, 0.06]
    widthOrientationRefSteps = numpy.asarray(widthOrientationRefSteps)
    widthOrientationRefSteps = widthOrientationRefSteps/180*numpy.pi
    widthCenterRefSteps = [0.10, 0.06, 0.04]
    distanceThreshold = 8
    distanceThresholds = [5, 4, 3]
    boxWidth = 48
    
    
    ### FOLDERS ###
    processingFolder = './Output_r%s/UnassembledImageProcessing'%runNumber        
    if not os.path.exists(processingFolder):
        os.mkdir(processingFolder)
    processingFiguresFolder = '%s/BgSubtractionAndPeakSearchPlots'%processingFolder
    if not os.path.exists(processingFiguresFolder):
        os.mkdir(processingFiguresFolder)
        
    ### EXTRACT GEOMETRY ###
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    #plotGeometry.plotGeometry(xGeometry_np, yGeometry_np)
    
    ### EXTRACT LATTICES ###    
    latticeObjectsPath = './Output_r%s/OrientationAndCellRefinement/r%s_refinedLatticesDictionary.pkl'%(runNumber, runNumber)
    fLattices = open(latticeObjectsPath, 'rb')
    latticesDictionary = pickle.load(fLattices)
    fLattices.close()
    
    ### LOOP ON LATTICES ###
    for myKey, myLattice in latticesDictionary.items():
        processFlag = 0
        if not imageSelection == '' and not latticeSelection == 0:
            if myLattice.runNumber == runNumber and myLattice.imageNumber == imageSelection and myLattice.latticeNumberInImage == latticeSelection:
                processFlag = 1
                bgPlanePlotFlag = 1
                singleSpotFigureFlag = 1
        else:
            if myLattice.runNumber == runNumber:
                processFlag = 1
                bgPlanePlotFlag = 0
                singleSpotFigureFlag = 0
            
        
        if processFlag == 1:
            print "Run %s Image %s Lattice %s"%(myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage)
            
            imageCenter = [0, 0]
            imageCenter = numpy.asarray(imageCenter, dtype=numpy.float32)    
            myLattice.imageCenter = imageCenter
            
            tiltAngle = float(myLattice.tiltAngle)/180*numpy.pi
            waveVector = 2*numpy.pi/myLattice.wavelength
            detectorDistance = myLattice.detectorDistance
            pixelSize = myLattice.pixelSize
            
            ### LISTS TO FOLLOW LATTICE REFINEMENT BEHAVIOUR ###
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
            fOpen = open('%s/%s-Processing_r%s_img%s_lattice%s.txt'%(processingFolder, minimizationMethod,
                                                                     myLattice.runNumber, 
                                                                     myLattice.imageNumber, 
                                                                     myLattice.latticeNumberInImage), 'w') 
            fOpen.write('Run: %s \nImage: %s \nLattice: %s\n'%(myLattice.runNumber, myLattice.imageNumber, myLattice.latticeNumberInImage))
            fOpen.write('START: \nLattice orientation: %.5f \nCell size: %.5f \n'%(myLattice.refinedInPlaneOrientation, myLattice.refinedCellSize)) 
            fOpen.write('Center x: %.6f\n'%myLattice.imageCenter[0])
            fOpen.write('Center y: %.6f\n'%myLattice.imageCenter[1])
            fOpen.write('\nConnected peaks detection:\n\n')  
            fOpen.write('    n    h    k    Predicted x    Predicted y    Local noise     Detected x     Detected y     Detected I    Distance (pxls)\n\n')
            
            ### LOOP ON PREDICTED PEAKS ###
            print 'Raw sectors extraction, background subtraction and connected peaks detection.\n'
            startTime = time.time()
            nPredictedPeak = 0         
            nDetected = 0
            detectedPeaks = []
            detectedPeaksDictionary = {}
            refinedPredictedPattern = myLattice.refinedPredictedPattern
            
            for predictedPeak in refinedPredictedPattern:
                if predictedPeak[4] <= lowResLimit and predictedPeak[4] >= highResLimit:
                    nPredictedPeak = nPredictedPeak + 1
                    
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
                        localNoise = bgSubtractionDictionary['localNoise']
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
                            moduleRotation = numpy.arctan(deltaYgeo/deltaXgeo) 
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
                        peakDetectionMask = peakDetectionDictionary['peakDetectionMask']
                        xCenterOfMass_vector = peakDetectionDictionary['xCenterOfMass_vector'] # In lab frame, in m and relative to detector center, to be compared to xPredicted
                        yCenterOfMass_vector = peakDetectionDictionary['yCenterOfMass_vector'] # In lab frame, in m and relative to detector center, to be compared to yPredicted
                        integratedIntensity_vector = peakDetectionDictionary['integratedIntensity_vector']
                        distanceFromPrediction_vector = peakDetectionDictionary['distanceFromPrediction_vector']
                        
                        detectedPeakObject = spotClass.diffractionSpot(nPredictedPeak, predictedPeak[0], predictedPeak[1],  
                                                                       boxMatrix, maskBoxMatrix, bgSubtractedBoxMatrix, peakDetectionMask,
                                                                       xLeft, yDown)
                        detectedPeakObject.setBoxIndices(i_good-yDown, j_good-xLeft)
                        
                        nConnectedPeaks = len(integratedIntensity_vector)
                        if nConnectedPeaks > 0:
                            minIdx = distanceFromPrediction_vector.index(min(distanceFromPrediction_vector))
                            if distanceFromPrediction_vector[minIdx] <= distanceThreshold: ### KEEP THIS EXPERIMENTAL PEAK (relaxed threshold) ###
                                nDetected = nDetected + 1
                                detectedPeak = []
                                xCoM = xCenterOfMass_vector[minIdx]
                                yCoM = yCenterOfMass_vector[minIdx]
                                detectedPeak.append(predictedPeak[0])                        # 0 ---> h
                                detectedPeak.append(predictedPeak[1])                        # 1 ---> k
                                detectedPeak.append(xCoM)                                    # 2 ---> x CoM in m, in lab frame (to be compared with Xpredicted)
                                detectedPeak.append(yCoM)                                    # 3 ---> y CoM in m, in lab frame (to be compared with Ypredicted)
                                detectedPeak.append(integratedIntensity_vector[minIdx])      # 4 ---> Integrated intensity
                                detectedPeak.append(distanceFromPrediction_vector[minIdx])   # 5 ---> Discrepancy (pxls) between prediction and observation (center of mass)
                                detectedPeaks.append(detectedPeak)
                                                         
                                detectedPeakObject.xCoM = xCoM
                                detectedPeakObject.yCoM = yCoM
                                detectedPeakObject.connectedI = integratedIntensity_vector[minIdx]
                                detectedPeakObject.setDistanceFromPrediction(distanceFromPrediction_vector[minIdx])
                                
                                fOpen.write('%5d%5d%5d%15.8f%15.8f%15.2f%15.8f%15.8f%15.2f%19.2f\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1], 
                                                                                                       Xpredicted, Ypredicted, localNoise, xCoM, yCoM,                                                                                                                    
                                                                                                       integratedIntensity_vector[minIdx], distanceFromPrediction_vector[minIdx])) 
                            else:
                                fOpen.write('%5d%5d%5d \tDetected peak(s) too far from prediction.\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1]))
                        else:
                            fOpen.write('%5d%5d%5d \tNo detected peaks.\n\n'%(nPredictedPeak, predictedPeak[0], predictedPeak[1]))
                            
                        detectedPeaksDictionary['Spot%d'%nPredictedPeak] = detectedPeakObject       
                                                         
            detectedPeaks = numpy.asarray(detectedPeaks, dtype=numpy.float32)
                           
            fOpen.write('Detected peaks (initial distance from prediction below %d pxls): %d'%(distanceThreshold, nDetected))
            
            endTime = time.time() - startTime
            print 'It took %.2f s\n'%endTime
            
            nDetectedThreshold = nPredictedPeak * fractionDetectedThreshold
            print 'N predicted: %d'%nPredictedPeak
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
                if minimizationMethod == '4Dbf':
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
                    print result.fun
                    ### RESULTS ###
                    refinedOrientation = result.x[1]/285
                    refinedCellSize = result.x[0]/4
                    refinedCenterX = result.x[2]/10
                    refinedCenterY = result.x[3]/10
                
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
                      
                fOpen.write('\nIteration: %d'%nIteration )
                fOpen.write('\nLattice orientation: %.5f\n'%myLattice.refinedInPlaneOrientation)
                fOpen.write('Cell size: %.5f\n'%myLattice.refinedCellSize)
                fOpen.write('Center x: %.6f\n'%myLattice.imageCenter[0])
                fOpen.write('Center y: %.6f'%myLattice.imageCenter[1])
                
                # RECALCULATE DISTANCES #
                detectedPeaks = recalculateDistance.recalculateDistanceFunction(imageCenter, predictedPattern, detectedPeaks, pixelSize)
                nIteration = nIteration + 1
                ### END 4D REFINEMENT ###
                                
            fOpen.close()
            refinementEnd = time.time() - refinementStart
            print 'Refinement took %.2f s\n'%refinementEnd  
            
            #myLattice.detectedPeaks = detectedPeaks  
            myLattice.setRefinedPattern(predictedPattern)
            
            ### PLOT REFINEMENT BEHAVIOUR ###    
            myLattice.refinedLatticeOrientations = refinedLatticeOrientations
            myLattice.refinedCellSizes = refinedCellSizes
            myLattice.refinedCenterXs = refinedCenterXs
            myLattice.refinedCenterYs = refinedCenterYs
            myLattice.nDetectedAndMatchedPeaks = nDetectedAndMatchedPeaks
            myLattice.avgLatticeErrors = avgLatticeErrors    
            myLattice.refinementBehaviourPlot(processingFolder, minimizationMethod)
            
            
            ### INTEGRATE AND PRODUCE SINGLE SPOT FIGURES AND GLOBAL PATTERN FIGURE ###
            integratedPeaks = []
            spotsDictionary = {}
            
            for mySpotKey, mySpot in detectedPeaksDictionary.items():            
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
                        integratedPeaks.append(integratedPeak)
                        spotsDictionary['Spot%d'%mySpot.n] = mySpot
                    else:                                                                                                       # refined predicted position falls out of the box
                        mySpot.setFinalBoxIndices(-1, -1)
                else:                                                                                                           # successFlag == 0 refined predicted position not found in geometry
                    mySpot.setFinalBoxIndices(-1, -1)
                
                
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
            for i in range(0, len(integratedPeaks)):
                for j in range(i+1, len(integratedPeaks)):
                    if n_data[i] > n_data[j]:
                        nTemp = n_data[i]
                        hTemp = h_data[i]
                        kTemp = k_data[i]
                        qRodTemp = qRod_data[i]
                        iTemp = I_data[i]
                        icTemp = Ic_data[i]
                        
                        n_data[i] = n_data[j]
                        h_data[i] = h_data[j]
                        k_data[i] = k_data[j]
                        qRod_data[i] = qRod_data[j]
                        I_data[i] = I_data[j]
                        Ic_data[i] = Ic_data[j]
                        
                        n_data[j] = nTemp
                        h_data[j] = hTemp
                        k_data[j] = kTemp
                        qRod_data[j] = qRodTemp
                        I_data[j] = iTemp  
                        Ic_data[j] = icTemp
                        
            orderedIntegratedIntensities = numpy.zeros((len(integratedPeaks), 4))
            for i in range(0, len(integratedPeaks)):
                    orderedIntegratedIntensities[i,0]=h_data[i]
                    orderedIntegratedIntensities[i,1]=k_data[i]
                    orderedIntegratedIntensities[i,2]=qRod_data[i]
                    orderedIntegratedIntensities[i,3]=Ic_data[i]
                    
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
            
if __name__ == "__main__":
    print "\n**** CALLING processing ****"
    processing(sys.argv[1:])