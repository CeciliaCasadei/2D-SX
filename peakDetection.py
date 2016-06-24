# -*- coding: utf-8 -*-
import numpy

def peakDetectionFunction(bgSubtractedBoxMatrix, 
                          localNoise, pixelSize,
                          moduleRotation, 
                          labX_Origin, labY_Origin,
                          maskBorder_left, maskBorder_down,
                          Xpredicted, Ypredicted):
                              
    xCenterOfMass_vector = []
    yCenterOfMass_vector = []
    integratedIntensity_vector = []
    distanceFromPrediction_vector = []
    
    peakDetectionMask = numpy.zeros(shape = bgSubtractedBoxMatrix.shape)        
    i,j = numpy.unravel_index(bgSubtractedBoxMatrix.argmax(), bgSubtractedBoxMatrix.shape)
    maxI = bgSubtractedBoxMatrix[i,j]
    while maxI >= 8*localNoise:
        xConnected = []                             ### NEW CONNECTED PEAK ###
        xConnected.append(j)
        yConnected = []
        yConnected.append(i)
        iConnected = []
        iConnected.append(maxI)
        peakDetectionMask[i,j]=1
        
        searchFlag = 1
        while searchFlag == 1:                      ### SEARCH ON NEIGHBOURING PIXELS ###
            # BUILD peakDetectionMask_extend       
            peakDetectionMask_extend = numpy.zeros(shape = bgSubtractedBoxMatrix.shape)   # 48x48
            peakDetectionMask_nRows = peakDetectionMask.shape[0]   
            peakDetectionMask_nCols = peakDetectionMask.shape[1]
            peakDetectionMaskIndices = numpy.argwhere(peakDetectionMask == 1)
            for peakDetectionMaskIndex in peakDetectionMaskIndices:
                rowIdx = peakDetectionMaskIndex[0]
                columnIdx = peakDetectionMaskIndex[1]
                peakDetectionMask_extend[rowIdx, columnIdx] = 1
                peakDetectionMask_extend[max([0,rowIdx-1]), columnIdx] = 1
                peakDetectionMask_extend[min([peakDetectionMask_nRows-1, rowIdx+1]), columnIdx] = 1
                peakDetectionMask_extend[rowIdx, max([0, columnIdx-1])] = 1
                peakDetectionMask_extend[rowIdx, min([columnIdx+1, peakDetectionMask_nCols-1])] = 1                      
            # BUILD peakDetectionMask_difference
            peakDetectionMask_difference = numpy.zeros(shape = bgSubtractedBoxMatrix.shape) 
            peakDetectionMask_difference = peakDetectionMask_extend - peakDetectionMask
            # LOOP ON ELEMENTS OF peakDetectionMask_difference WITH VALUE OF 1
            pDM_differenceIndices = numpy.argwhere(peakDetectionMask_difference == 1)
            for pDM_differenceIndex in pDM_differenceIndices:
                rowIdx = pDM_differenceIndex[0]
                columnIdx = pDM_differenceIndex[1]
                intensityValue = bgSubtractedBoxMatrix[rowIdx, columnIdx]
                if intensityValue >= max([8*localNoise, maxI/1000]):
                    peakDetectionMask[rowIdx, columnIdx] = 1
                    xConnected.append(columnIdx)
                    yConnected.append(rowIdx)
                    iConnected.append(intensityValue)
                else:
                    peakDetectionMask_difference[rowIdx, columnIdx] = 0                            
            # IF peakDetectionMask_difference IS 0 EVERYWHERE, STOP SEARCH
            stopCheckIndices = numpy.argwhere(peakDetectionMask_difference == 1)
            if len(stopCheckIndices) == 0:
                searchFlag = 0
                integratedIntensity = numpy.sum(iConnected)
                integratedIntensity_vector.append(integratedIntensity)
                
                xI_Product = numpy.multiply(xConnected, iConnected).sum()
                xCoM = xI_Product / integratedIntensity             # float in (0, 96) or smaller.
                                
                yI_Product = numpy.multiply(yConnected, iConnected).sum()
                yCoM = yI_Product / integratedIntensity             # float in (0, 96) or smaller.            
                
                # TRANSLATION AND ROTATION FOR CONVERSION TO LAB FRAME #
                labX = labX_Origin + pixelSize*((xCoM-maskBorder_left)*numpy.cos(moduleRotation) - (yCoM-maskBorder_down)*numpy.sin(moduleRotation))
                labY = labY_Origin + pixelSize*((xCoM-maskBorder_left)*numpy.sin(moduleRotation) + (yCoM-maskBorder_down)*numpy.cos(moduleRotation))
                xCenterOfMass_vector.append(labX)
                yCenterOfMass_vector.append(labY)
                
                distanceFromPrediction = numpy.sqrt((labX-Xpredicted)**2 + (labY-Ypredicted)**2)/pixelSize # pxls
                distanceFromPrediction_vector.append(distanceFromPrediction)
                
                peakDetectionMask_inverse = 1 - peakDetectionMask
                peakDetectionBox_updated = numpy.multiply(bgSubtractedBoxMatrix, peakDetectionMask_inverse)
                i,j = numpy.unravel_index(peakDetectionBox_updated.argmax(), peakDetectionBox_updated.shape)
                maxI = bgSubtractedBoxMatrix[i,j]
                
    peakDetectionDictionary = {'xCenterOfMass_vector': xCenterOfMass_vector, 
                               'yCenterOfMass_vector': yCenterOfMass_vector,
                               'integratedIntensity_vector': integratedIntensity_vector, 
                               'distanceFromPrediction_vector': distanceFromPrediction_vector,
                               'peakDetectionMask': peakDetectionMask}
    return peakDetectionDictionary