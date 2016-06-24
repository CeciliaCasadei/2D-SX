# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 13:00:16 2016
@author: casadei_c
MEMBER OF diffractionImage CLASS
CALCULATE 256 TRIAL PREDICTED PATTERNS
CORRESPONDING TO 256 TRIAL ORIENTATIONS IN [0, 2pi]
STORE PREDICTED PATTERNS IN DICTIONARY: self.referencePredictedPattern
[h k qx qy dMin q azimuth rotatedAzimuth detectorAzimuth diffractionAngle detectorRadius]
"""
import pickle
import math
import numpy
import os

import calculatePredictedPattern

def getPredictedPatternFunction(self, detectorDistance, pixelSize, cellSize, trialInPlaneRotations):
    tiltAngle = float(self.tiltAngle)/180*numpy.pi
    waveVector = 2*math.pi/self.wavelength
    
    reciprocalLatticeFile = './Output_r%s/ReferenceReciprocalLattice/reciprocalLattice_cellSize_%.3f.pkl'%(self.runNumber, cellSize)
    
    if not os.path.exists(reciprocalLatticeFile):
        print 'File %s not found.'%reciprocalLatticeFile
        
    else:
        fRead = open(reciprocalLatticeFile, 'rb')
        reciprocalLattice = pickle.load(fRead)
        fRead.close()
        
        predictedPatternsDictionary = {}
        n = 0        
        for trialInPlaneRotation in trialInPlaneRotations:                     # 256 trial in-plane orientations
            n = n + 1
            predictedPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, waveVector, tiltAngle, detectorDistance, pixelSize)
            predictedPatternsDictionary['%s'%n] = predictedPattern             # Keys: 1 -> 256
        self.referencePredictedPattern = predictedPatternsDictionary