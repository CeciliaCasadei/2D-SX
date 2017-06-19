# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:18:03 2016
@author: casadei_c
LOCAL BACKGROUND SUBTRACTION
LOW FLUCTUATION POINTS ARE SELECTED AND THEIR INTENSITY IS USED TO MODEL THE BACKGROUND BY
- COMPUTING THE AVERAGE
OR
- FITTING WITH A PLANE
"""
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import scipy.linalg

import plotBgPlane

def backGroundSubtractionFunction(boxMatrix, maskBoxMatrix, bgSubtractionMethod, 
                                  fileName, runNumber, imageNumber, latticeNumberInImage, 
                                  nPeak, plotFlag, directoryName = ''):                                   
    matplotlib.pyplot.ioff
    
    gridStep = 6
    subBoxWidth = 6
    lowFluctuationThreshold = 2 # 2.5 in MATLAB
    
    xGrid = range(6, boxMatrix.shape[1], gridStep)   # 6 12 18 ... 90 (max 15 values)  IT WAS range(0, boxMatrix.shape[1]+6, gridStep)
    yGrid = range(6, boxMatrix.shape[0], gridStep)   # 6 12 18 ... 90 (max 15 values)  IT WAS range(0, boxMatrix.shape[1]+6, gridStep)
    stdDevMatrix = numpy.zeros((len(yGrid), len(xGrid))) # 15x15 or smaller
    avgsMatrix = numpy.zeros((len(yGrid), len(xGrid)))   # 15x15 or smaller
    stdDevVector = []
    boxMatrixNrows = boxMatrix.shape[0]
    boxMatrixNcolumns = boxMatrix.shape[1]
    xGridItemIdx = 0        
    for xGridItem in xGrid:
        yGridItemIdx = 0            
        for yGridItem in yGrid:
            stdDevMatrix[yGridItemIdx, xGridItemIdx] = 1000
            avgsMatrix[yGridItemIdx, xGridItemIdx] = -1
            
            xLeft_subBox = max([xGridItem - subBoxWidth, 0])
            xRight_subBox = min([xGridItem + subBoxWidth, boxMatrixNcolumns])
            yDown_subBox = max(yGridItem - subBoxWidth, 0)
            yUp_subBox = min(yGridItem + subBoxWidth, boxMatrixNrows)
            
            subBoxMatrix = boxMatrix[yDown_subBox:yUp_subBox, xLeft_subBox:xRight_subBox]          # 12x12 or smaller
            subMaskBoxMatrix = maskBoxMatrix[yDown_subBox:yUp_subBox, xLeft_subBox:xRight_subBox]  # 12x12 or smaller
            myValidIndices = numpy.argwhere(subMaskBoxMatrix == 1)                                 # Exclude neighbouring modules
            nValid = myValidIndices.shape[0]
            
            if nValid > 1 and nValid > (subMaskBoxMatrix.shape[0] * subMaskBoxMatrix.shape[1] / 2):                                        
                avg = numpy.sum(subBoxMatrix) / nValid              
                mySquare = numpy.square(subBoxMatrix)
                mySum = numpy.sum(mySquare)                                                        # print mySum.dtype ---> float64
                myVariance = (mySum / nValid) - (avg ** 2)
                if myVariance > 0:
                    stdDev = numpy.sqrt(myVariance)
                    stdDevVector.append(stdDev)
                    stdDevMatrix[yGridItemIdx, xGridItemIdx] = stdDev
                    avgsMatrix[yGridItemIdx, xGridItemIdx] = avg                    
                else:
                    print "BG SUBTRACTION PROBLEM: VARIANCE %.18f"%myVariance                    
            yGridItemIdx = yGridItemIdx + 1
        xGridItemIdx = xGridItemIdx + 1
    localNoise = numpy.percentile(stdDevVector, 15)
    
    xSample = []
    ySample = []
    intensitySample = []           
    
    lowFluctuationIndices = numpy.argwhere(stdDevMatrix <= lowFluctuationThreshold * localNoise)
    for myIndices in lowFluctuationIndices:
        myRowIdx = myIndices[0]
        myColumnIdx = myIndices[1]
        xSample.append(xGrid[myColumnIdx])
        ySample.append(yGrid[myRowIdx])
        intensitySample.append(avgsMatrix[myRowIdx, myColumnIdx])
    
    if bgSubtractionMethod == 'average':
        ### METHOD AVERAGE ###
        backGround = numpy.average(intensitySample)        
    else:
        ### METHOD BACKGROUND PLANE ###
        myData = numpy.c_[xSample, ySample, intensitySample]       
        A = numpy.c_[myData[:,0], myData[:,1], numpy.ones(myData.shape[0])]
        C,_,_,_ = scipy.linalg.lstsq(A, myData[:,2]) 
        
        myX, myY = numpy.meshgrid(numpy.arange(0, boxMatrix.shape[1], 1), numpy.arange(0, boxMatrix.shape[0], 1))
        backGround = C[0]*myX + C[1]*myY + C[2]          
        if plotFlag == 1:
            plotBgPlane.plotBgPlane(myX, myY, backGround,
                                    myData,
                                    C,
                                    directoryName, runNumber, imageNumber, latticeNumberInImage, nPeak)
    ### BG SUBTRACTION ###
    bgSubtractedBoxMatrix = boxMatrix - backGround
                                                          
    return {'localNoise':localNoise, 'bgSubtractedBoxMatrix':bgSubtractedBoxMatrix}