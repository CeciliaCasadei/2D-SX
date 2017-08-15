# -*- coding: utf-8 -*-
import numpy
import scipy
import scipy.interpolate
import scipy.optimize
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
    

    
def integrate(spotMatrix, expansion_factor=0):
    
    # PARAMETERS
    integrationRadius = 5*(10**expansion_factor)
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
    

    
def integrate_ellipse(spotMatrix, sigma_x, sigma_y, multiplicative_factor, expansion_factor=0):
    
    # PARAMETERS
    nCountsPerPhoton = 26
    
    # N COUNTS TO N PHOTONS CONVERSION
    spotMatrix =  spotMatrix/nCountsPerPhoton
    
    # PREPARE INTEGRATION MASK (ELLIPSE)
    integrationMask = numpy.zeros((spotMatrix.shape))   
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    x_axis = multiplicative_factor * sigma_x * (10**expansion_factor)
    y_axis = multiplicative_factor * sigma_y * (10**expansion_factor)
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
    

def poly_4(x, offset, a, b):
    y = offset + a*x**2 + b*x**4
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
    
def calculateBackground_nBg(imageSector, lowFluctuationThreshold=2):
        
    # PARAMETERS
    gridStep = 3
    subBoxWidth = 3
    
    # CALCULATE STD DEVIATION AND MEAN INTENSITY ON A GRID, CALCULATE SECTOR LOCAL NOISE
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

    sectorSize = imageSector.shape[0] * imageSector.shape[1]
    n_gridPoints = len(xGrid) * len(yGrid)
    n_fitPoints = len(intensitySample)
    if n_fitPoints > n_gridPoints:
        print 'PROBLEM!!!'
    n_bgPixels = int(sectorSize * ( float(n_fitPoints) / n_gridPoints))      
                               
    return backGround, n_bgPixels
    
    

def calculate_detectorAzimuth(xGeometry_np, yGeometry_np, i, j):   
    xDetector = xGeometry_np[i, j]
    yDetector = yGeometry_np[i, j]
    if xDetector != 0:
        detectorAzimuth = numpy.arctan(yDetector/xDetector)
        if xDetector < 0 and yDetector > 0:
            detectorAzimuth = detectorAzimuth + numpy.pi
        if xDetector < 0 and yDetector < 0:
            detectorAzimuth = detectorAzimuth - numpy.pi
    else:
        if yDetector > 0:
            detectorAzimuth = + numpy.pi /2
        else:
            detectorAzimuth = - numpy.pi /2
    return detectorAzimuth


    
def calculate_moduleRotation(xGeometry_np, yGeometry_np, i, j):
    deltaXgeo = float(xGeometry_np[i, j+1] - xGeometry_np[i, j])
    deltaYgeo = float(yGeometry_np[i, j+1] - yGeometry_np[i, j])

    if deltaXgeo != 0:
        moduleRotation = numpy.arctan(deltaYgeo/deltaXgeo) # From module to lab frame
        if deltaXgeo < 0 and deltaYgeo > 0:
            moduleRotation = moduleRotation + numpy.pi
        elif deltaXgeo < 0 and deltaYgeo < 0:
            moduleRotation = moduleRotation - numpy.pi
    else:
        if deltaYgeo > 0:
            moduleRotation = numpy.pi / 2
        else:
            moduleRotation = -numpy.pi / 2
    return moduleRotation


    
def clockWiseRotation(spotMatrix, rotationAngle):
    x_windows_pix = range(-spotMatrix.shape[1]/2, +spotMatrix.shape[1]/2)  # -18, -17, ..., +17
    y_windows_pix = range(-spotMatrix.shape[0]/2, +spotMatrix.shape[0]/2)  # -18, -17, ..., +17

    [X_windows_pix, Y_windows_pix] = numpy.meshgrid(x_windows_pix, y_windows_pix)
    X_windows_pix = numpy.asarray(X_windows_pix, dtype=numpy.float32)
    Y_windows_pix = numpy.asarray(Y_windows_pix, dtype=numpy.float32)
    
    X_windows_pix_rotated = numpy.cos(rotationAngle)*X_windows_pix - numpy.sin(rotationAngle)*Y_windows_pix
    Y_windows_pix_rotated = numpy.sin(rotationAngle)*X_windows_pix + numpy.cos(rotationAngle)*Y_windows_pix
    
    f = scipy.interpolate.interp2d(x_windows_pix, y_windows_pix, spotMatrix, kind='linear')
    spotMatrix_rotated = numpy.zeros(spotMatrix.shape)
    
    for columnIndex in range(0, spotMatrix_rotated.shape[1]):
        for rowIndex in range(0, spotMatrix_rotated.shape[0]):
            rotated_x = X_windows_pix_rotated[rowIndex, columnIndex]
            rotated_y = Y_windows_pix_rotated[rowIndex, columnIndex]
            rotated_f = f(rotated_x, rotated_y)
            spotMatrix_rotated[rowIndex, columnIndex] = rotated_f   
            
    return spotMatrix_rotated



def do_gaussFit(sector):
    try:                    
        n_x = sector.shape[1]
        n_y = sector.shape[0]
        x = numpy.linspace(0, n_x-1, n_x)
        y = numpy.linspace(0, n_y-1, n_y)
        x, y = numpy.meshgrid(x, y)                      # column_index, row_index
        
        data = sector.ravel()           
        data = data.T
        data = numpy.asarray(data)
        data = data.flatten()
        
        initial_x = float(n_x)/2
        initial_y = float(n_y)/2
        initial_guess = (numpy.amax(sector), initial_x, initial_y, 1.0, 1.0)
        
        popt, pcov = scipy.optimize.curve_fit(twoD_Gaussian_simple, (x, y), data, p0=initial_guess)
        data_fitted = twoD_Gaussian_simple((x, y), *popt)       
        
        refined_amplitude = popt[0]
        refined_x0 = popt[1]
        refined_y0 = popt[2]
        refined_sigma_x = popt[3]
        refined_sigma_y = popt[4]
              
        ### ANALYTICAL GAUSSIAN INTEGRAL ###
        gauss_integral = 2 * numpy.pi * refined_amplitude * refined_sigma_x * refined_sigma_y
        
    except:
        print 'Gaussian fit not possible'
        gauss_integral = numpy.nan
        refined_amplitude = numpy.nan
        refined_x0 = numpy.nan
        refined_y0 = numpy.nan
        refined_sigma_x = numpy.nan
        refined_sigma_y = numpy.nan
        data = numpy.nan
        data_fitted = numpy.nan
        
    return refined_sigma_x, refined_sigma_y, refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted
    
    
    
def do_gaussFit_fixed_sigmas(sector, sigX, sigY):
    try:                    
        n_x = sector.shape[1]
        n_y = sector.shape[0]
        x = numpy.linspace(0, n_x-1, n_x)
        y = numpy.linspace(0, n_y-1, n_y)
        x, y = numpy.meshgrid(x, y)                      # column_index, row_index
        
        data = sector.ravel()           
        data = data.T
        data = numpy.asarray(data)
        data = data.flatten()
        
        initial_x = float(n_x)/2
        initial_y = float(n_y)/2
        initial_guess = (numpy.amax(sector), initial_x, initial_y)
                      
        popt, pcov = scipy.optimize.curve_fit(lambda (x, y), amplitude, xo, yo: twoD_Gaussian_simple((x, y), amplitude, xo, yo, sigX, sigY), (x, y), data, p0=initial_guess)   
        
        refined_amplitude = popt[0]
        refined_x0 = popt[1]
        refined_y0 = popt[2]
        
        data_fitted = twoD_Gaussian_simple((x, y), refined_amplitude, refined_x0, refined_y0, sigX, sigY)
              
        ### ANALYTICAL GAUSSIAN INTEGRAL ###
        gauss_integral = 2 * numpy.pi * refined_amplitude * sigX * sigY
        
    except:
        print 'Gaussian fit not possible'
        gauss_integral = numpy.nan
        refined_amplitude = numpy.nan
        refined_x0 = numpy.nan
        refined_y0 = numpy.nan
        data = numpy.nan
        data_fitted = numpy.nan
        
    return refined_x0, refined_y0, refined_amplitude, gauss_integral, data, data_fitted



def recenter(sector, x0, y0, precision_factor, truncated_halfWidth):
    
    precision_factor = 10**precision_factor                      # 10
    truncated_halfWidth = precision_factor*truncated_halfWidth   # 150
    
    expanded_sector = numpy.zeros((precision_factor*sector.shape[0], precision_factor*sector.shape[1])) # 500x500
    x0 = (precision_factor*x0)+(precision_factor/2) # NOT: x0 = precision_factor*x0
    y0 = (precision_factor*y0)+(precision_factor/2) # NOT: y0 = precision_factor*y0
    x0 = int(round(x0))
    y0 = int(round(y0))
    for i in range(0, sector.shape[0]):
        for j in range(0, sector.shape[1]):
            element = sector[i, j]
            for m in range(i*precision_factor, i*precision_factor+precision_factor):
                for n in range(j*precision_factor, j*precision_factor+precision_factor):
                    expanded_sector[m, n] = element
    recentered_sum = expanded_sector[y0-truncated_halfWidth:y0+truncated_halfWidth, x0-truncated_halfWidth:x0+truncated_halfWidth] 
    return recentered_sum # 300x300
    
    
    
def translate(spotMatrix, dX, dY):
    # Positive dX is a translation to the left

    x_windows_pix = range(0, +spotMatrix.shape[1])  # -18, -17, ..., +17
    y_windows_pix = range(0, +spotMatrix.shape[0])  # -18, -17, ..., +17
    
    [X_windows_pix, Y_windows_pix] = numpy.meshgrid(x_windows_pix, y_windows_pix) # [columnIdx, rowIdx]
    X_windows_pix = numpy.asarray(X_windows_pix, dtype=numpy.float32)
    Y_windows_pix = numpy.asarray(Y_windows_pix, dtype=numpy.float32)
    
    X_windows_pix_translated = X_windows_pix+dX
    Y_windows_pix_translated = Y_windows_pix+dY
   
    f = scipy.interpolate.interp2d(x_windows_pix, y_windows_pix, spotMatrix, kind='linear')
    spotMatrix_translated = numpy.zeros(spotMatrix.shape)
    
    for columnIndex in range(0, spotMatrix_translated.shape[1]):
        for rowIndex in range(0, spotMatrix_translated.shape[0]):
            translated_x = X_windows_pix_translated[rowIndex, columnIndex]
            translated_y = Y_windows_pix_translated[rowIndex, columnIndex]
            translated_f = f(translated_x, translated_y)
            spotMatrix_translated[rowIndex, columnIndex] = translated_f
            
    return spotMatrix_translated
    
    

def buildMasks(sector, precision_factor, multiplicative_factor, sigma_x, sigma_y):
    # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
    integrationMask = numpy.zeros((sector.shape))   
    ringMask        = numpy.zeros((sector.shape)) 
    
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    x_axis = 10**precision_factor * multiplicative_factor * sigma_x
    y_axis = 10**precision_factor * multiplicative_factor * sigma_y
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
    
    integrationMask[numpy.where(distance < 1)] = 1
    ringMask[numpy.where(distance > 1.5)] = 1
    ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
    
    return integrationMask, ringMask