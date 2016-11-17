# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 15:07:13 2016
@author: casadei_c
CALCULATE PREDICTED DIFFRACTION PATTERN 
(GIVEN LATTICE ORIENTATION AND CELL SIZE, THE LATTER DETERMINES THE RECIPROCAL LATTICE INPUT).
FORMULAS VALID FOR GENERAL TILT.
"""
import numpy
def calculatePredictedPatternFunction(reciprocalLattice, trialInPlaneRotation, waveVector, tiltAngle, detectorDistance, pixelSize):
    predictedPattern = numpy.zeros((len(reciprocalLattice), 13))
        
    index = 0
    for i in reciprocalLattice:                                                # Reciprocal lattice vectors
        h = i[0]
        k = i[1]
        qx = i[2]
        qy = i[3]
        resolution = i[4]
        q = i[5]
        
        azimuth = numpy.arcsin(qy/q)                                           # Between -pi/2 and +pi/2
        if qx < 0:
            azimuth = numpy.pi - azimuth                                       # Between pi/2 and 3pi/2
        if azimuth < 0:
            azimuth = 2*numpy.pi + azimuth                                     # Azimuth always positive
        
        rotatedAzimuth = azimuth + trialInPlaneRotation
        if rotatedAzimuth > 0:
            rotatedAzimuth = rotatedAzimuth % (2*numpy.pi)
        else:
            rotatedAzimuth = - rotatedAzimuth
            rotatedAzimuth = rotatedAzimuth % (2*numpy.pi)
            rotatedAzimuth = - rotatedAzimuth
            rotatedAzimuth = 2*numpy.pi + rotatedAzimuth                       # RotatedAzimuth in [0, 2pi]
            
        qrod = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q/waveVector)**2 - 2*q/waveVector*numpy.sin(tiltAngle)*numpy.sin(rotatedAzimuth)))
        qxy = numpy.sqrt(q**2 + numpy.sin(tiltAngle)**2 * (qrod**2 - q**2 * numpy.sin(rotatedAzimuth)**2) + 2*q*qrod * numpy.sin(rotatedAzimuth) * numpy.cos(tiltAngle) * numpy.sin(tiltAngle) )
        sinDetectorAzimuth = (q*numpy.sin(rotatedAzimuth)*numpy.cos(tiltAngle) + qrod*numpy.sin(tiltAngle)) /qxy
                
        if abs(sinDetectorAzimuth  - 1) < 0.000001:
            sinDetectorAzimuth = 1.0
        if abs(sinDetectorAzimuth  + 1) < 0.000001:
            sinDetectorAzimuth = -1.0                
                
        detectorAzimuth = numpy.arcsin(sinDetectorAzimuth)                     # Between -pi/2 and +pi/2
              
        if rotatedAzimuth > numpy.pi/2 and rotatedAzimuth < 3*numpy.pi/2:
            detectorAzimuth = numpy.pi - detectorAzimuth                       # Between -pi/2 and +3pi/2
                
        if detectorAzimuth < 0:
            detectorAzimuth = 2*numpy.pi + detectorAzimuth                     # Between 0 and 2pi
                
        if sinDetectorAzimuth > 1 or sinDetectorAzimuth < -1:
            print 'Problem: h=%d k=%d qx=%.2f qy = %.2f q=%.2f rotated azimuth=%.2f qxy=%.2f sin of detector azimuth=%.15f  detector azimuth=%.3f'%(h,k,qx,qy,q,rotatedAzimuth,qxy,sinDetectorAzimuth,detectorAzimuth)
               
        diffractionAngle = 2*numpy.arcsin((numpy.sqrt(qrod**2+q**2))/(2*waveVector)) # Between -pi and +pi
        detectorRadius = detectorDistance * numpy.tan(diffractionAngle) / pixelSize  # Pxls
        
        ### LP CORRECTION ###
        Pfactor = 1 - (numpy.sin(diffractionAngle))**2 * (numpy.cos(detectorAzimuth))**2
        xDetector = pixelSize * detectorRadius * numpy.cos(detectorAzimuth)    # m
        yDetector = pixelSize * detectorRadius * numpy.sin(detectorAzimuth)    # m
        N = numpy.sqrt(xDetector**2 + yDetector**2 + detectorDistance**2)
        Lfactor = numpy.multiply([0, -numpy.sin(tiltAngle), numpy.cos(tiltAngle)], [xDetector/N, yDetector/N, detectorDistance/N])
        Lfactor = Lfactor.sum()
        LPfactor = Lfactor / Pfactor
           
        predictedPattern[index, 0] = h                                         # int
        predictedPattern[index, 1] = k                                         # int
        predictedPattern[index, 2] = qx                                        # A^(-1)
        predictedPattern[index, 3] = qy                                        # A^(-1)
        predictedPattern[index, 4] = resolution                                # A
        predictedPattern[index, 5] = q                                         # A^(-1)
        predictedPattern[index, 6] = azimuth                                   # radians
        predictedPattern[index, 7] = rotatedAzimuth                            # radians
        predictedPattern[index, 8] = detectorAzimuth                           # radians, [0, 2pi]
        predictedPattern[index, 9] = diffractionAngle                          # radians, [-pi, pi]
        predictedPattern[index, 10] = detectorRadius                           # pxls
        predictedPattern[index, 11] = qrod                                     # A^(-1)
        predictedPattern[index, 12] = LPfactor                                 # float
    
        index = index + 1
    
    return predictedPattern