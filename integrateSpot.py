# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:24:35 2016
Member of diffractionSpot class
@author: casadei_c
"""
import numpy

def integrate(self):
    
    bgSubtractedBoxMatrixPhotons =  self.bgSubtractedBoxMatrix/self.nCountsPerPhoton
    
    integrationMask = numpy.zeros((bgSubtractedBoxMatrixPhotons.shape))
   
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    radialDistance = numpy.sqrt((rowIdx-self.iFinal)**2+(colIdx-self.jFinal)**2)
    integrationMask[numpy.where(radialDistance < self.integrationRadius)] = 1
    self.integrationMask = integrationMask
    
    integrationFlag = 0
    if 0 <= self.iFinal - self.integrationRadius  and self.iFinal + self.integrationRadius < self.bgSubtractedBoxMatrix.shape[0]:
        if 0 <= self.jFinal - self.integrationRadius  and self.jFinal + self.integrationRadius < self.bgSubtractedBoxMatrix.shape[1]:
            if numpy.multiply(self.maskBoxMatrix, integrationMask).sum() == integrationMask.sum():
                integrationFlag = 1
    
    self.integrationFlag = integrationFlag        
    
    if integrationFlag == 1:
        integratedIntensity = numpy.multiply(integrationMask, bgSubtractedBoxMatrixPhotons).sum()
    else:
        integratedIntensity = numpy.nan
    return integratedIntensity