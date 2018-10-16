# -*- coding: utf-8 -*-
import numpy

def buildMasks(sector, multiplicative_factor, sigma_x, sigma_y):
    # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
    integrationMask = numpy.zeros((sector.shape))   
    ringMask        = numpy.zeros((sector.shape)) 
    
    colIdx, \
    rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), 
                            numpy.arange(integrationMask.shape[0]))
    x_axis = multiplicative_factor * sigma_x
    y_axis = multiplicative_factor * sigma_y
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
    
    integrationMask[numpy.where(distance < 1)] = 1
    ringMask[numpy.where(distance > 1.5)] = 1
    ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
    
    return integrationMask, ringMask
