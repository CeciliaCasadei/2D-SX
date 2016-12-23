# -*- coding: utf-8 -*-
import numpy
import scipy
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D



def twoD_Gaussian_simple((x, y), amplitude, xo, yo, sigma_x, sigma_y):    
    g = amplitude * numpy.exp( - ((x-xo)**2)/(2*sigma_x**2) - ((y-yo)**2)/(2*sigma_y**2)   )
    return g.ravel()
    

    
def integrate(spotMatrix):
    
    # PARAMETERS
    integrationRadius = 5
    nCountsPerPhoton = 26
    
    # N COUNTS TO N PHOTONS CONVERSION
    spotMatrix =  spotMatrix/nCountsPerPhoton
    
    # PREPARE INTEGRATION MASK (CIRCLE)
    integrationMask = numpy.zeros((spotMatrix.shape))   
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    radialDistance = numpy.sqrt((rowIdx-centerY)**2 + (colIdx-centerX)**2)
    integrationMask[numpy.where(radialDistance < integrationRadius)] = 1
    
    # INTEGRATE BY SUMMING ON CIRCLE
    integratedIntensity = numpy.multiply(integrationMask, spotMatrix).sum()
    
    return integratedIntensity
    

    
def integrate_ellipse(spotMatrix, sigma_x, sigma_y, multiplicative_factor):
    
    # PARAMETERS
    nCountsPerPhoton = 26
    #multiplicative_factor = 4.0 #2.5
    
    # N COUNTS TO N PHOTONS CONVERSION
    spotMatrix =  spotMatrix/nCountsPerPhoton
    
    # PREPARE INTEGRATION MASK (ELLIPSE)
    integrationMask = numpy.zeros((spotMatrix.shape))   
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    x_axis = multiplicative_factor * sigma_x
    y_axis = multiplicative_factor * sigma_y
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
    integrationMask[numpy.where(distance < 1)] = 1
    
    # INTEGRATE BY SUMMING ON ELLIPSE
    integratedIntensity = numpy.multiply(integrationMask, spotMatrix).sum()
    
    return integratedIntensity

    

def line(x, offset, slope):
    y = offset + slope*x 
    return y


    
def line_plus_sigmoid(x, offset, slope, x0, k, scale):
    y = offset + slope*x + scale / (1 + numpy.exp(-k*(x-x0)))
    return y
    

     
def quadratic(x, offset, a, b):
    y = offset + a*x + b*x**2
    return y


def quadratic_no_first_order(x, offset, a):
    y = offset + a*x**2
    return y
    

def quadratic_plus_sigmoid(x, offset, a, x0, k, scale):
    y = offset + a*x**2 + scale / (1 + numpy.exp(-k*(x-x0)))
    return y

    
def Correlate(x1, x2):
    x1Avg = numpy.average(x1)
    x2Avg = numpy.average(x2)
    numTerm = numpy.multiply(x1-x1Avg, x2-x2Avg)
    num = numTerm.sum()
    resX1Sq = numpy.multiply(x1-x1Avg, x1-x1Avg)
    resX2Sq = numpy.multiply(x2-x2Avg, x2-x2Avg)
    den = numpy.sqrt(numpy.multiply(resX1Sq.sum(), resX2Sq.sum()))
    CC = num/den
    return CC    



def calculateBackground(imageSector, folder, label, h, k):
    
    # PARAMETERS
    gridStep = 3
    subBoxWidth = 3
    lowFluctuationThreshold = 2 # 2.5 in MATLAB
    
    # CALCULATE STD DEVIATIN AND MEAN INTENSITY ON A GRID, CALCULATE SECTOR LOCAL NOISE
    xGrid = range(3, imageSector.shape[1], gridStep)     # 3 6 9 ... 27 (9 values) 
    yGrid = range(3, imageSector.shape[0], gridStep)     # 3 6 9 ... 27 (9 values) 
    stdDevMatrix = numpy.zeros((len(yGrid), len(xGrid))) # 9x9
    avgsMatrix   = numpy.zeros((len(yGrid), len(xGrid)))   # 9x9
    stdDevVector = []
    sectorNrows = imageSector.shape[0]
    sectorNcolumns = imageSector.shape[1]
    xGridItemIdx = 0        
    for xGridItem in xGrid:
        yGridItemIdx = 0            
        for yGridItem in yGrid:
            
            xLeft_subBox = max([xGridItem - subBoxWidth, 0])
            xRight_subBox = min([xGridItem + subBoxWidth, sectorNcolumns])
            yDown_subBox = max(yGridItem - subBoxWidth, 0)
            yUp_subBox = min(yGridItem + subBoxWidth, sectorNrows)
            
            subBoxMatrix = imageSector[yDown_subBox:yUp_subBox, xLeft_subBox:xRight_subBox] # 6x6
                                         
            avg = numpy.sum(subBoxMatrix) / (subBoxMatrix.shape[0] * subBoxMatrix.shape[1])        
            mySquare = numpy.square(subBoxMatrix)
            mySum = numpy.sum(mySquare)                                                      
            myVariance = (mySum / (subBoxMatrix.shape[0] * subBoxMatrix.shape[1])  ) - (avg ** 2)
            if myVariance > 0:
                stdDev = numpy.sqrt(myVariance)
                stdDevVector.append(stdDev)
                stdDevMatrix[yGridItemIdx, xGridItemIdx] = stdDev
                avgsMatrix[yGridItemIdx, xGridItemIdx] = avg                    
            else:
                print "*************BG SUBTRACTION PROBLEM: VARIANCE %.18f**************"%myVariance                    
            yGridItemIdx = yGridItemIdx + 1
        xGridItemIdx = xGridItemIdx + 1
    localNoise = numpy.percentile(stdDevVector, 15)
    
    # EXTRACT INTENSITIES OF LOW FLUCTUATION POINTS FOR SUBSEQUENT BG PLANE CALCULATION
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

    # BACKGROUND PLANE FIT
    myData = numpy.c_[xSample, ySample, intensitySample]       
    A = numpy.c_[myData[:,0], myData[:,1], numpy.ones(myData.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, myData[:,2]) 
    
    myX, myY = numpy.meshgrid(numpy.arange(0, imageSector.shape[1], 1), numpy.arange(0, imageSector.shape[0], 1))
    backGround = C[0]*myX + C[1]*myY + C[2]        
    
    # PLOT BACKGROUND PLANE
    if not os.path.exists('%s/BackgroundPlane'%folder):
        os.mkdir('%s/BackgroundPlane'%folder)
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(myData[:,0], myData[:,1], myData[:,2], c='r', marker='o', s=6)
    mySurface = myAxes.plot_surface(myX, myY, backGround, rstride=1, cstride=1, cmap=cm.coolwarm,
                                    linewidth=0, antialiased=False, alpha = 0.3)
    myFigure.colorbar(mySurface, shrink=0.5, aspect=5)
  
    if C[1] >= 0:
        s1 = '+'
    else:
        s1 = '-'
    if C[2] >= 0:
        s2 = '+'
    else:
        s2 = '-'        
    matplotlib.pyplot.title('z = %.3fx %s %.3fy %s %.2f'%(C[0], s1, abs(C[1]), s2, abs(C[2])), fontsize = 10, y=1.08)
    
    myAxes.zaxis.set_major_locator(LinearLocator(10))
    myAxes.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    myAxes.set_zlim(0, 40)
    myAxes.set_xlabel('X')
    myAxes.set_ylabel('Y')
    myAxes.set_zlabel('Z')  
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
   
    matplotlib.pyplot.savefig('%s/BackgroundPlane/BgPlane_%s_orbit_%d_%d.png'
                               %(folder, label, h, k))   
    matplotlib.pyplot.close()                           
    return backGround
    
    
    
def calculateBackground_noImg(imageSector):
    
    # PARAMETERS
    gridStep = 3
    subBoxWidth = 3
    lowFluctuationThreshold = 2 # 2.5 in MATLAB
    
    # CALCULATE STD DEVIATIN AND MEAN INTENSITY ON A GRID, CALCULATE SECTOR LOCAL NOISE
    xGrid = range(3, imageSector.shape[1], gridStep)     # 3 6 9 ... 27 (9 values) 
    yGrid = range(3, imageSector.shape[0], gridStep)     # 3 6 9 ... 27 (9 values) 
    stdDevMatrix = numpy.zeros((len(yGrid), len(xGrid))) # 9x9
    avgsMatrix   = numpy.zeros((len(yGrid), len(xGrid)))   # 9x9
    stdDevVector = []
    sectorNrows = imageSector.shape[0]
    sectorNcolumns = imageSector.shape[1]
    xGridItemIdx = 0        
    for xGridItem in xGrid:
        yGridItemIdx = 0            
        for yGridItem in yGrid:
            
            xLeft_subBox = max([xGridItem - subBoxWidth, 0])
            xRight_subBox = min([xGridItem + subBoxWidth, sectorNcolumns])
            yDown_subBox = max(yGridItem - subBoxWidth, 0)
            yUp_subBox = min(yGridItem + subBoxWidth, sectorNrows)
            
            subBoxMatrix = imageSector[yDown_subBox:yUp_subBox, xLeft_subBox:xRight_subBox] # 6x6
                                         
            avg = numpy.sum(subBoxMatrix) / (subBoxMatrix.shape[0] * subBoxMatrix.shape[1])        
            mySquare = numpy.square(subBoxMatrix)
            mySum = numpy.sum(mySquare)                                                      
            myVariance = (mySum / (subBoxMatrix.shape[0] * subBoxMatrix.shape[1])  ) - (avg ** 2)
            if myVariance > 0:
                stdDev = numpy.sqrt(myVariance)
                stdDevVector.append(stdDev)
                stdDevMatrix[yGridItemIdx, xGridItemIdx] = stdDev
                avgsMatrix[yGridItemIdx, xGridItemIdx] = avg                    
            else:
                print "*************BG SUBTRACTION PROBLEM: VARIANCE %.18f**************"%myVariance                    
            yGridItemIdx = yGridItemIdx + 1
        xGridItemIdx = xGridItemIdx + 1
    localNoise = numpy.percentile(stdDevVector, 15)
    
    # EXTRACT INTENSITIES OF LOW FLUCTUATION POINTS FOR SUBSEQUENT BG PLANE CALCULATION
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

    # BACKGROUND PLANE FIT
    myData = numpy.c_[xSample, ySample, intensitySample]       
    A = numpy.c_[myData[:,0], myData[:,1], numpy.ones(myData.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, myData[:,2]) 
    
    myX, myY = numpy.meshgrid(numpy.arange(0, imageSector.shape[1], 1), numpy.arange(0, imageSector.shape[0], 1))
    backGround = C[0]*myX + C[1]*myY + C[2]         
                               
    return backGround
