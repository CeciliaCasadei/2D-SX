# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 15:20:55 2016
@author: casadei_c
MEMBER OF diffractionImage CLASS
IDENTIFY AND INDEX MULTIPLE LATTICES IN ONE IMAGE
OBJECTS OF THE CLASS Lattice ARE GENERATED
"""
import os
import numpy
import pickle

import latticeClass

def indexingFunction(self, detectorDistance, pixelSize, 
                     radialTolerance, pixelTolerance, azimuthTolerance, 
                     minNofPeaksPerLattice, 
                     referenceCellSize,
                     geometryFile):   
                         
    nInPlaneAngles = 256
    trialInPlaneAngles = numpy.arange(-numpy.pi, numpy.pi, 2*numpy.pi/nInPlaneAngles)              # 256 elements
    self.getPredictedPattern(detectorDistance, pixelSize, referenceCellSize, trialInPlaneAngles)   # Calculate 256 trial predicted patterns 
                                                                                                   # and store result in self.referencePredictedPattern
                                                                                                   # (dictionary with 256 items)
                                                                                                   # [h k qx qy dMin q azimuth rotatedAzimuth detectorAzimuth 
                                                                                                   # diffractionAngle detectorRadius qRod LPfactor]
                                                                                                   # NB: azimuth on detector in [0, 2pi]    
    nValidPeaks = 0
    for i in range(0, self.nPeaks):
        if self.orderedPeaksMatrix[i,5] == 1:
            nValidPeaks = nValidPeaks + 1
       
    nPredictedSpots = 0
    for i in self.referencePredictedPattern['1']:
        nPredictedSpots = nPredictedSpots + 1                                  # 432
    
    latticeDictionaryFromOneImage = {}
    nLatticesInImage = 0
    startIdx = 0
    while nValidPeaks > minNofPeaksPerLattice and startIdx < self.nPeaks:
        
        # Look for most intense, valid experimental peak.
        for i in range(startIdx, self.nPeaks):
            startIdx = startIdx +1
            if self.orderedPeaksMatrix[i,5] == 1:
                intenseExpPeakIdx = i
                break
        
        # Extract radius and azimuth of most intense valid experimental peak
        intenseExpPeakRadius  = self.orderedPeaksMatrix[intenseExpPeakIdx, 3]  # Radius, pxls
        intenseExpPeakAzimuth = self.orderedPeaksMatrix[intenseExpPeakIdx, 4]  # Between 0 and 2pi
        
        # For every trial orientation, for every predicted spot radially close to intense peak, 
        # report value of azimuth difference to intense peak in deltaAzimuthMatrix.
        deltaAzimuthMatrix = numpy.zeros((nPredictedSpots, nInPlaneAngles))    # 432 x 256 
        for i in range(0,nPredictedSpots):
            for j in range(0, nInPlaneAngles):
                deltaAzimuthMatrix[i,j] = 500
        
        for inPlaneAngle in range(1, nInPlaneAngles+1):                                        # 1 to 256 included
            predictionInOneOrientation = self.referencePredictedPattern['%s'%inPlaneAngle]     # 432 x 12
            for predictedSpot in range(0, nPredictedSpots):                                    # 0 to 431 included              
                predictedRadius  = predictionInOneOrientation[predictedSpot,10]                # Radius, pxls
                predictedAzimuth = predictionInOneOrientation[predictedSpot,8]                 # Between 0 and 2pi
                if numpy.abs(predictedRadius-intenseExpPeakRadius) <= radialTolerance:   
                    phiTolerance = min([float(azimuthTolerance)/180*numpy.pi, float(pixelTolerance)/predictedRadius])
                    if numpy.abs(predictedAzimuth-intenseExpPeakAzimuth) <= phiTolerance:
                        deltaAzimuthMatrix[predictedSpot, inPlaneAngle-1] = numpy.abs(predictedAzimuth - intenseExpPeakAzimuth)
                    elif predictedAzimuth <= phiTolerance/2 and intenseExpPeakAzimuth >= 2*numpy.pi - phiTolerance/2:
                        deltaAzimuthMatrix[predictedSpot, inPlaneAngle-1] = predictedAzimuth + 2*numpy.pi - intenseExpPeakAzimuth
                    elif intenseExpPeakAzimuth <= phiTolerance/2 and predictedAzimuth >= 2*numpy.pi - phiTolerance/2:
                        deltaAzimuthMatrix[predictedSpot, inPlaneAngle-1] = intenseExpPeakAzimuth + 2*numpy.pi - predictedAzimuth
                    else:
                        continue
        
        # Find minimum of deltaAzimuthMatrix -> possibleRotationAngle            
        i,j = numpy.unravel_index(deltaAzimuthMatrix.argmin(),deltaAzimuthMatrix.shape)
        possibleRotationAngleIndex = j                                                          # From 0
        possibleRotationAngle = trialInPlaneAngles[j]
        
        # Count n of experimental peaks matching predicted peaks in the selected possibleRotationAngle
        nMatches = 0
        myKey = possibleRotationAngleIndex + 1
        possiblePredictedPattern = self.referencePredictedPattern['%s'%myKey]
        for i in range(0, self.nPeaks):
            if self.orderedPeaksMatrix[i,5] == 1:
                expRadius  = self.orderedPeaksMatrix[i,3]                      # Radius, pxls
                expAzimuth = self.orderedPeaksMatrix[i,4]                      # Azimuth in [0, 2pi]
                for j in range(0, nPredictedSpots):
                    predictedRadius  = possiblePredictedPattern[j,10]          # Radius, pxls
                    predictedAzimuth = possiblePredictedPattern[j,8]           # Azimuth in [0, 2pi]
                    if numpy.abs(predictedRadius-expRadius) <= radialTolerance:
                        phiTolerance = min([float(azimuthTolerance)/180*numpy.pi, float(pixelTolerance)/predictedRadius])
                        if numpy.abs(predictedAzimuth-expAzimuth) <= phiTolerance:
                            nMatches = nMatches + 1
                        elif predictedAzimuth <= phiTolerance/2 and expAzimuth >= 2*numpy.pi - phiTolerance/2:
                            nMatches = nMatches + 1
                        elif expAzimuth <= phiTolerance/2 and predictedAzimuth >= 2*numpy.pi - phiTolerance/2:
                            nMatches = nMatches + 1
                    else:
                        continue
        
        # If the number of matches is big enough, possibleRotationangle defines a Lattice orientation
        if nMatches >= minNofPeaksPerLattice:
            nLatticesInImage = nLatticesInImage + 1            
            rotationAngle = possibleRotationAngle
            
            # Store all matches in indexedPeaksTable
            indexedPeaksTable = numpy.zeros((nMatches, 11))
            indexedPeaksTableRow = 0
            for i in range(0, self.nPeaks):
                if self.orderedPeaksMatrix[i,5] == 1:
                    expRadius  = self.orderedPeaksMatrix[i,3]                  # Radius, pxls
                    expAzimuth = self.orderedPeaksMatrix[i,4]                  # Azimuth in [0, 2pi]
                    for j in range(0, nPredictedSpots):
                        predictedRadius  = possiblePredictedPattern[j,10]
                        predictedAzimuth = possiblePredictedPattern[j,8]
                        if numpy.abs(predictedRadius-expRadius) <= radialTolerance and indexedPeaksTableRow < nMatches:
                            phiTolerance = min([float(azimuthTolerance)/180*numpy.pi, float(pixelTolerance)/predictedRadius])
                            indexedPeaksTable[indexedPeaksTableRow,0] = possiblePredictedPattern[j,0]                     # h
                            indexedPeaksTable[indexedPeaksTableRow,1] = possiblePredictedPattern[j,1]                     # k
                            indexedPeaksTable[indexedPeaksTableRow,2] = expRadius                                         # experimental radius
                            indexedPeaksTable[indexedPeaksTableRow,3] = expAzimuth                                        # experimental azimuth
                            indexedPeaksTable[indexedPeaksTableRow,4] = self.orderedPeaksMatrix[i,2]                      # experimental intensity
                            indexedPeaksTable[indexedPeaksTableRow,5] = abs(expRadius - predictedRadius)                  # radial difference
                            indexedPeaksTable[indexedPeaksTableRow,7] = i                                                 # experimental peak n                                               
                            indexedPeaksTable[indexedPeaksTableRow,8] = j                                                 # predicted peak n
                            indexedPeaksTable[indexedPeaksTableRow,9]  = predictedRadius                                  # predicted radius
                            indexedPeaksTable[indexedPeaksTableRow,10] = predictedAzimuth                                 # predicted azimuth
                            if numpy.abs(predictedAzimuth-expAzimuth) <= phiTolerance:                              
                                self.orderedPeaksMatrix[i,5] = 0
                                indexedPeaksTable[indexedPeaksTableRow,6] = abs(expAzimuth - predictedAzimuth)            # azimuth difference
                                indexedPeaksTableRow = indexedPeaksTableRow + 1
                            elif predictedAzimuth <= phiTolerance/2 and expAzimuth >= 2*numpy.pi - phiTolerance/2:
                                self.orderedPeaksMatrix[i,5] = 0
                                indexedPeaksTable[indexedPeaksTableRow,6] = predictedAzimuth + 2*numpy.pi - expAzimuth
                                indexedPeaksTableRow = indexedPeaksTableRow + 1
                            elif expAzimuth <= phiTolerance/2 and predictedAzimuth >= 2*numpy.pi - phiTolerance/2:    
                                self.orderedPeaksMatrix[i,5] = 0
                                indexedPeaksTable[indexedPeaksTableRow,6] = expAzimuth + 2*numpy.pi - predictedAzimuth
                                indexedPeaksTableRow = indexedPeaksTableRow + 1
                            else:
                                continue
                                
            # GENERATE Lattice OBJECT, STORE IN DICTIONARY 
            latticeObject = latticeClass.Lattice(self.fileName, self.imageNumber, self.runNumber, self.tiltAngle, 
                                                 rotationAngle, self.wavelength, pixelSize, detectorDistance, nMatches, 
                                                 referenceCellSize, indexedPeaksTable, nLatticesInImage)
            latticeDictionaryFromOneImage['%s_Lattice_%d'%(self.fileName, nLatticesInImage)] = latticeObject           
            
            # Update n of valid experimental peaks.
            nValidPeaks = 0
            for i in range(0, self.nPeaks):
                if self.orderedPeaksMatrix[i,5] == 1:
                    nValidPeaks = nValidPeaks + 1
    
    # Count n of lattices in the image.
    nLattices = 0
    for i,j in latticeDictionaryFromOneImage.items():
        nLattices = nLattices + 1
    
    # Save individual image lattices file.    
    if nLattices == 0:
        print 'Image %s - %s:\tno lattices found.'%(self.imageNumber.zfill(5), self.fileName)
    else:    
        if not os.path.exists('./Output_r%s/LatticeIndexing'%self.runNumber):   
            os.mkdir('./Output_r%s/LatticeIndexing'%self.runNumber)          
        outputFile = './Output_r%s/LatticeIndexing/latticeDictionary_r%s_image_%s.pkl'%(self.runNumber, self.runNumber, self.imageNumber)  
        f = open(outputFile, 'wb')
        pickle.dump(latticeDictionaryFromOneImage, f)
        f.close()