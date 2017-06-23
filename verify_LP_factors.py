# -*- coding: utf-8 -*-
import numpy

import buildReciprocalLattice
import calculatePredictedPattern

cellSize = 62.45
hmax = 100
kmax = 100
resolutionLimit = 4

reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(cellSize, hmax, kmax, resolutionLimit)

trialInPlaneRotation = 0
wavelength = 1.48
wavevector = 2*numpy.pi / wavelength
tiltAngle = 0
detectorDistance = 0.235     # m
pixelSize = 0.000110         # m

predictedPattern = calculatePredictedPattern.calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, wavevector, tiltAngle, detectorDistance, pixelSize)

fSummary = open('./predictedPattern.txt', 'w')
fSummary.write('  h    k        q_x        q_y      d_min        q   azimuth    rotated azimuth    azimuth on detector    diffraction angle    radius on detector      qRod      LP\n')
for predictedPatternLine in predictedPattern:
    fSummary.write('%3d%5d%11.4f%11.4f%11.2f%9.3f%10.3f%19.4f%23.4f%21.4f%22.4f%10.5f%8.4f\n'
                    %(predictedPatternLine[0], predictedPatternLine[1], 
                      predictedPatternLine[2], predictedPatternLine[3], 
                      predictedPatternLine[4], predictedPatternLine[5], 
                      predictedPatternLine[6], predictedPatternLine[7], 
                      predictedPatternLine[8], predictedPatternLine[9], 
                      predictedPatternLine[10], predictedPatternLine[11], 
                      predictedPatternLine[12]))
fSummary.close() 