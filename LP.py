# -*- coding: utf-8 -*-
import numpy

def get_LP(xDetector, 
           yDetector,
           q_2D,
           qRod,
           waveVector,
           detectorDistance,
           tiltAngle):
    
    cosDetectorAzimuth = xDetector/numpy.sqrt(xDetector**2 + yDetector**2)     
    diffractionAngle = 2*numpy.arcsin((numpy.sqrt(qRod**2+q_2D**2))/(2*waveVector)) # Between -pi and +pi
    Pfactor = 1 - (numpy.sin(diffractionAngle))**2 * cosDetectorAzimuth**2
    N = numpy.sqrt(xDetector**2 + yDetector**2 + detectorDistance**2)
    Lfactor = numpy.multiply([0,           -numpy.sin(tiltAngle), numpy.cos(tiltAngle)], 
                             [xDetector/N, yDetector/N,           detectorDistance/N  ])
    Lfactor = Lfactor.sum()
    LPfactor = Lfactor / Pfactor
    
    return LPfactor
