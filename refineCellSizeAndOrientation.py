# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 17:30:43 2016
@author: casadei_c
refineCellSizeAndOrientation
MEMBER OF Lattice CLASS
CALLED BY orientationAndCellRefinement
GIVEN A LIST OF CELL SIZE REFINEMENT STEPS AND A LIST OF IN-PLANE ORIENTATION REFINEMENT STEPS,
THE LATTICE ERROR IS CALCULATED IN EACH GRID POINT OF THE 2D REFINEMENT PARAMETERS SPACE AND MINIMIZED.
"""
import os
import numpy
import pickle

import buildReciprocalLattice

import calculatePredictedPattern
import calculateLatticeErrorForRefinement     
      
def refineCellSizeAndOrientation(self, 
                                 sizeRefinementSteps, 
                                 orientationRefinementSteps, 
                                 folderName, 
                                 resolutionLimit):
    
    orientationRefinementSteps = self.inPlaneRotation + orientationRefinementSteps

    tiltAngle = float(self.tiltAngle)/180*numpy.pi
    waveVector = 2*numpy.pi/self.wavelength
      
    latticeErrorMatrix = numpy.zeros((len(orientationRefinementSteps), 
                                      len(sizeRefinementSteps)))
    nSizeStep = 0
    for sizeRefinementStep in sizeRefinementSteps: # (21) trial unit cell sizes
        nSizeStep = nSizeStep + 1     
        
#        reciprocalLatticeDataFile = '%s/reciprocalLattice_cellSize_%.3f.pkl'%(RLfolderName, sizeRefinementStep)
#        if os.path.exists(reciprocalLatticeDataFile):
#            fRead = open(reciprocalLatticeDataFile, 'rb')
#            reciprocalLattice = pickle.load(fRead)
#            fRead.close()
#        else:
#            print 'File %s not found!'%reciprocalLatticeDataFile
#            return
        reciprocalLattice = \
        buildReciprocalLattice.buildReciprocalLatticeFunction(sizeRefinementStep, 
                                                              100, 
                                                              100, 
                                                              resolutionLimit)
        
        nRotationStep = 0       
        for trialInPlaneRotation in orientationRefinementSteps: # (21) trial in-plane orientations
            nRotationStep = nRotationStep + 1
            myPredictedPattern = \
            calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, 
                                                                        trialInPlaneRotation, 
                                                                        waveVector, 
                                                                        tiltAngle, 
                                                                        self.detectorDistance, 
                                                                        self.pixelSize)
            myPredictedPattern = numpy.asarray(myPredictedPattern, dtype=numpy.float32)
            indexedPeaksTable = numpy.asarray(self.indexedPeaksTable, dtype=numpy.float32)
            
            latticeError, \
            nMatchedPeaks = \
            calculateLatticeErrorForRefinement.calculateLatticeError(myPredictedPattern, 
                                                                     indexedPeaksTable)
                                                                     
            latticeErrorMatrix[nRotationStep-1,nSizeStep-1] = latticeError
            
    i,j = numpy.unravel_index(latticeErrorMatrix.argmin(),latticeErrorMatrix.shape)
    minError = latticeErrorMatrix[i, j]
    avgMinError = minError / nMatchedPeaks
    refinedOrientation = orientationRefinementSteps[i]
    refinedCellSize = sizeRefinementSteps[j]
    
    refinedRL = \
    buildReciprocalLattice.buildReciprocalLatticeFunction(refinedCellSize, 
                                                          100, 
                                                          100, 
                                                          resolutionLimit)
                                                              
#    reciprocalLatticeFile = '%s/reciprocalLattice_cellSize_%.3f.pkl'%(RLfolderName, refinedCellSize)
#    fRead = open(reciprocalLatticeFile, 'rb')
#    refinedRL = pickle.load(fRead)
#    fRead.close()
    
    refinedPattern = \
    calculatePredictedPattern.calculatePredictedPatternFunction(refinedRL, 
                                                                refinedOrientation, 
                                                                waveVector, 
                                                                tiltAngle, 
                                                                self.detectorDistance, 
                                                                self.pixelSize)
    
    return {'refinedOrientation': refinedOrientation, 
            'refinedCellSize': refinedCellSize, 
            'latticeErrorMatrix': latticeErrorMatrix, 
            'refinedPredictedPattern': refinedPattern,
            'averageError': avgMinError}