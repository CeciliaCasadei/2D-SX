# -*- coding: utf-8 -*-
import numpy
import warnings
import matplotlib.pyplot
import math
import os

def buildReciprocalLatticeFunction(cellSize, transformationMatrix, directVectors, hmax, kmax, resolutionLimit):
    directCell = cellSize * directVectors * transformationMatrix
    print directCell
    reciprocalCellRows = 2* numpy.pi * directCell.I
    reciprocalLattice = []
    i = 0
    for h in range(-hmax, +hmax+1):
        for k in range(-kmax, +kmax+1):
            reciprocalVector = [h, k]*reciprocalCellRows
            q_x = reciprocalVector[0,0]
            q_y = reciprocalVector[0,1]
            q = numpy.sqrt(q_x**2 + q_y**2)
            if q != 0:
                resolution = 2* numpy.pi / q
            if q != 0 and resolution >= resolutionLimit:
                reciprocalLatticeItem = []
                reciprocalLatticeItem.append(int(h)) 
                reciprocalLatticeItem.append(int(k))
                reciprocalLatticeItem.append(q_x)
                reciprocalLatticeItem.append(q_y)
                reciprocalLatticeItem.append(resolution)
                reciprocalLatticeItem.append(q)
                reciprocalLattice.append(reciprocalLatticeItem)
                i = i + 1
    
    # SORTING REFLECTIONS IN ORDER OF ASCENDING q:
    for i_sort in range(0, i):
        for j_sort in range(i_sort+1, i):
            if reciprocalLattice[i_sort][5] > reciprocalLattice[j_sort][5]:
                h = reciprocalLattice[i_sort][0]
                k = reciprocalLattice[i_sort][1]
                qx = reciprocalLattice[i_sort][2]
                qy = reciprocalLattice[i_sort][3]
                resolution = reciprocalLattice[i_sort][4]
                q = reciprocalLattice[i_sort][5]
                reciprocalLattice[i_sort][0] = reciprocalLattice[j_sort][0]
                reciprocalLattice[i_sort][1] = reciprocalLattice[j_sort][1]
                reciprocalLattice[i_sort][2] = reciprocalLattice[j_sort][2]
                reciprocalLattice[i_sort][3] = reciprocalLattice[j_sort][3]
                reciprocalLattice[i_sort][4] = reciprocalLattice[j_sort][4]
                reciprocalLattice[i_sort][5] = reciprocalLattice[j_sort][5]
                reciprocalLattice[j_sort][0] = h
                reciprocalLattice[j_sort][1] = k
                reciprocalLattice[j_sort][2] = qx
                reciprocalLattice[j_sort][3] = qy
                reciprocalLattice[j_sort][4] = resolution
                reciprocalLattice[j_sort][5] = q
    # END SORTING    
    return reciprocalLattice
    
def calculatePredictedPatternFunction(RLdata, trialInPlaneRotation, waveVector, tiltAngle, detectorDistance, pixelSize):
    predictedPattern = numpy.zeros((len(RLdata), 13))
        
    index = 0
    for i in RLdata: # Reciprocal lattice vectors
        h = i[0]
        k = i[1]
        qx = i[2]
        qy = i[3]
        resolution = i[4]
        q = i[5]
        azimuth = numpy.arcsin(qy/q)        #between -pi/2 and +pi/2
        if qx < 0:
            azimuth = numpy.pi - azimuth    #between pi/2 and 3pi/2
        if azimuth < 0:
            azimuth = 2*numpy.pi + azimuth
        # azimuth always positive
        rotatedAzimuth = azimuth + trialInPlaneRotation
        if rotatedAzimuth > 0:
            rotatedAzimuth = rotatedAzimuth % (2*numpy.pi)
        else:
            rotatedAzimuth = - rotatedAzimuth
            rotatedAzimuth = rotatedAzimuth % (2*numpy.pi)
            rotatedAzimuth = - rotatedAzimuth
            rotatedAzimuth = 2*numpy.pi + rotatedAzimuth
        # rotatedAzimuth in [0, 2pi]    
        qrod = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q/waveVector)**2 - 2*q/waveVector*numpy.sin(tiltAngle)*numpy.sin(rotatedAzimuth)))
        qxy = numpy.sqrt(q**2 + numpy.sin(tiltAngle)**2 * (qrod**2 - q**2 * numpy.sin(rotatedAzimuth)**2) + 2*q*qrod * numpy.sin(rotatedAzimuth) * numpy.cos(tiltAngle) * numpy.sin(tiltAngle) )
        sinDetectorAzimuth = (q*numpy.sin(rotatedAzimuth)*numpy.cos(tiltAngle) + qrod*numpy.sin(tiltAngle)) /qxy
                
        if abs(sinDetectorAzimuth  - 1) < 0.000001:
            sinDetectorAzimuth = 1.0

        if abs(sinDetectorAzimuth  + 1) < 0.000001:
            sinDetectorAzimuth = -1.0                
                
        detectorAzimuth = numpy.arcsin(sinDetectorAzimuth)          #between -pi/2 and +pi/2
              
        if rotatedAzimuth > numpy.pi/2 and rotatedAzimuth < 3*numpy.pi/2:
            detectorAzimuth = numpy.pi - detectorAzimuth            #between -pi/2 and +3pi/2
                
        if detectorAzimuth < 0:
            detectorAzimuth = 2*numpy.pi + detectorAzimuth          #between 0 and 2pi
                
        if sinDetectorAzimuth > 1 or sinDetectorAzimuth < -1:
            print 'Problem: h=%d k=%d qx=%.2f qy = %.2f q=%.2f rotated azimuth=%.2f qxy=%.2f sin of detector azimuth=%.15f  detector azimuth=%.3f'%(h,k,qx,qy,q,rotatedAzimuth,qxy,sinDetectorAzimuth,detectorAzimuth)
               
        diffractionAngle = 2*numpy.arcsin((numpy.sqrt(qrod**2+q**2))/(2*waveVector)) #between -pi and +pi
        detectorRadius = detectorDistance * numpy.tan(diffractionAngle) / pixelSize
           
        predictedPattern[index, 0] = h
        predictedPattern[index, 1] = k
        predictedPattern[index, 2] = qx
        predictedPattern[index, 3] = qy
        predictedPattern[index, 4] = resolution
        predictedPattern[index, 5] = q
        predictedPattern[index, 6] = azimuth
        predictedPattern[index, 7] = rotatedAzimuth
        predictedPattern[index, 8] = detectorAzimuth
        predictedPattern[index, 9] = diffractionAngle
        predictedPattern[index, 10] = detectorRadius
        predictedPattern[index, 11] = qrod
        
        ### LP CORRECTION ###
        Pfactor = 1 - (numpy.sin(diffractionAngle))**2 * (numpy.cos(detectorAzimuth))**2
        xDetector = pixelSize * detectorRadius * numpy.cos(detectorAzimuth)     # m
        yDetector = pixelSize * detectorRadius * numpy.sin(detectorAzimuth)     # m
        N = numpy.sqrt(xDetector**2 + yDetector**2 + detectorDistance**2)
        Lfactor = numpy.multiply([0, -numpy.sin(tiltAngle), numpy.cos(tiltAngle)], [xDetector/N, yDetector/N, detectorDistance/N])
        Lfactor = Lfactor.sum()
        LPfactor = Lfactor / Pfactor
        print Lfactor
        print Pfactor
        print LPfactor
        print '\n'
        predictedPattern[index, 12] = LPfactor
        
        index = index + 1
    
    return predictedPattern    

def logPredictedPattern(predictedPattern, myString):
    if not os.path.exists('./patternSymmetry'):
        os.mkdir('./patternSymmetry')   
    fLog = open('./patternSymmetry/%s.log'%myString, 'w')
    fLog.write('    h     k  xDetector  yDetector         qRod\n')
    for diffractionSpot in predictedPattern:
        detectorX = diffractionSpot[10] * math.cos(diffractionSpot[8])
        detectorY = diffractionSpot[10] * math.sin(diffractionSpot[8])
        fLog.write('%5d %5d %10.3f %10.3f %12.5f\n'%(diffractionSpot[0], diffractionSpot[1], detectorX, detectorY, diffractionSpot[11]))
    fLog.close()
    
def plotPredictedPattern(predictedPattern1, predictedPattern2, resolutionRadii, detectorDistance, wavelength, myString):
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.close()
    myDPI = 96 

    myFigure = matplotlib.pyplot.figure(figsize=(25, 25), dpi=myDPI, facecolor = 'w')   # Figure object
    myAxes = myFigure.add_subplot(1,1,1) # Axes object
    
#    myX1 = []
#    myY1 = []
#    myLabels1 = []    
#    for myRow1 in predictedPattern1:
#        h_idx1 = myRow1[0]
#        h_idx1 = int(h_idx1)
#        k_idx1 = myRow1[1]
#        k_idx1 = int(k_idx1)
#        hk1 = '%d, %d'%(h_idx1, k_idx1)
#        predictedRadius1 = myRow1[10]
#        predictedAzimuth1 = myRow1[8]
#        predicted_x1 = predictedRadius1 * math.cos(predictedAzimuth1)
#        predicted_y1 = predictedRadius1 * math.sin(predictedAzimuth1)
#        
#        myX1.append(predicted_x1)
#        myY1.append(predicted_y1)
#        myLabels1.append(hk1)
#        
#    myAxes.scatter(myX1, myY1, color='b', marker="o", facecolors='none', s=70)
#    for label1, x1, y1 in zip(myLabels1, myX1, myY1):
#        matplotlib.pyplot.annotate(label1, xy = (x1, y1), xytext = (-3, 3), size = 18, 
#                                   textcoords = 'offset points', ha = 'right', va = 'bottom', color='b')
                                   
    
    
    myX2 = []
    myY2 = []
    myLabels2 = []    
    for myRow2 in predictedPattern2:
        h_idx2 = myRow2[0]
        h_idx2 = int(h_idx2)
        k_idx2 = myRow2[1]
        k_idx2 = int(k_idx2)
        hk2 = '%d, %d'%(h_idx2, k_idx2)
        predictedRadius2 = myRow2[10]
        predictedAzimuth2 = myRow2[8]
        predicted_x2 = predictedRadius2 * math.cos(predictedAzimuth2)
        predicted_y2 = predictedRadius2 * math.sin(predictedAzimuth2)
        
        myX2.append(predicted_x2)
        myY2.append(predicted_y2)
        myLabels2.append(hk2)
        
    myAxes.scatter(myX2, myY2, color='r', marker="o", facecolors='none', s=120)
    for label2, x2, y2 in zip(myLabels2, myX2, myY2):
        matplotlib.pyplot.annotate(label2, xy = (x2, y2), xytext = (+3, -3), size = 30, 
                                   textcoords = 'offset points', ha = 'left', va = 'top', color='b')
                                   
                                                                  
    myRadii = []
    for i in resolutionRadii:
        i = float(i)
        myRadius = detectorDistance/pixelSize*math.tan(2*math.asin(wavelength/(2*i)))
        circle = matplotlib.pyplot.Circle((0,0),myRadius,linewidth=0.5, color='b',fill=False)
        myAxes.add_artist(circle)
        myRadii.append(myRadius)
    
    matplotlib.pyplot.axhline(y=0, xmin=-870, xmax=870, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-870, ymax=870, linewidth=0.5, color = 'b')
    myAxes.set_title("%s \nBlue: identity, Red: transformed \nResolution circles: %.1f A, %.1f A, %.1f A"%(myString, resolutionRadii[0], resolutionRadii[1], resolutionRadii[2]), 
                     y=1.02, fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    
    myAxes.set_xlim([-max(myRadii)-20,+max(myRadii)+20])
    myAxes.set_ylim([-max(myRadii)-20,+max(myRadii)+20])
    
    myAxes.set_xlabel("Detector pxls, x", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("Detector pxls, y", fontsize = 22, rotation = 'vertical')
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    
    if not os.path.exists('./patternSymmetry'):
        os.mkdir('./patternSymmetry')   
    myFigure.savefig("./patternSymmetry/lowRes_%s.png"%myString)
    matplotlib.pyplot.close()


cellSize = 62.45
resolutionLimit = 20
hmax = 100
kmax = 100

directVectors = numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]])
directVectorsP3 = numpy.matrix([[-0.5, -numpy.sqrt(3)/2], [numpy.sqrt(3)/2, -0.5]]) * directVectors # rotated direct vectors
directVectorsP3P3 = numpy.matrix([[-0.5, numpy.sqrt(3)/2], [-numpy.sqrt(3)/2, -0.5]]) * directVectors # rotated direct vectors

identityMatrix = numpy.matrix([[1, 0],[0, 1]])
permutationMatrix = numpy.matrix([[0, 1],[1, 0]])
inversionMatrix = numpy.matrix([[-1, 0],[0, -1]])

reciprocalLattice_I = buildReciprocalLatticeFunction(cellSize, identityMatrix, directVectors, hmax, kmax, resolutionLimit)
reciprocalLattice_P = buildReciprocalLatticeFunction(cellSize, permutationMatrix, directVectors, hmax, kmax, resolutionLimit)
reciprocalLattice_i = buildReciprocalLatticeFunction(cellSize, inversionMatrix, directVectors, hmax, kmax, resolutionLimit)
reciprocalLattice_P3 = buildReciprocalLatticeFunction(cellSize, identityMatrix, directVectorsP3, hmax, kmax, resolutionLimit)
reciprocalLattice_P3P3 = buildReciprocalLatticeFunction(cellSize, identityMatrix, directVectorsP3P3, hmax, kmax, resolutionLimit)

inPlaneRotation = 0 * numpy.pi / 180
waveLength = 1.486567
waveVector = 2*numpy.pi / waveLength
tiltAngle = 0 * numpy.pi / 180
detectorDistance = 0.235
pixelSize = 0.000110
resolutionRadii = [50, 10, 7]    
predictedPattern_I = calculatePredictedPatternFunction(reciprocalLattice_I, 
                                                       inPlaneRotation, 
                                                       waveVector, 
                                                       tiltAngle, 
                                                       detectorDistance, 
                                                       pixelSize)
predictedPattern_P = calculatePredictedPatternFunction(reciprocalLattice_P,
                                                       inPlaneRotation, 
                                                       waveVector, 
                                                       tiltAngle, 
                                                       detectorDistance, 
                                                       pixelSize)
predictedPattern_i = calculatePredictedPatternFunction(reciprocalLattice_i,
                                                       inPlaneRotation, 
                                                       waveVector, 
                                                       tiltAngle, 
                                                       detectorDistance, 
                                                       pixelSize)
predictedPattern_P3 = calculatePredictedPatternFunction(reciprocalLattice_P3, 
                                                        inPlaneRotation, 
                                                        waveVector, 
                                                        tiltAngle, 
                                                        detectorDistance, 
                                                        pixelSize)
predictedPattern_P3P3 = calculatePredictedPatternFunction(reciprocalLattice_P3P3,
                                                          inPlaneRotation, 
                                                          waveVector, 
                                                          tiltAngle, 
                                                          detectorDistance, 
                                                          pixelSize)

plotPredictedPattern(predictedPattern_I, predictedPattern_I, resolutionRadii, detectorDistance, waveLength, 'DirectAxesIdentity')
#logPredictedPattern(predictedPattern_I, 'Identity')
plotPredictedPattern(predictedPattern_I, predictedPattern_P, resolutionRadii, detectorDistance, waveLength, 'DirectAxesPermutation')
#logPredictedPattern(predictedPattern_P, 'IndicesPermutation')
plotPredictedPattern(predictedPattern_I, predictedPattern_i, resolutionRadii, detectorDistance, waveLength, 'DirectAxesInversion')
#logPredictedPattern(predictedPattern_i, 'IndicesInversion')
plotPredictedPattern(predictedPattern_I, predictedPattern_P3, resolutionRadii, detectorDistance, waveLength, 'DirectAxesP3rotation')
#logPredictedPattern(predictedPattern_P3, 'P3rotation')
plotPredictedPattern(predictedPattern_I, predictedPattern_P3P3, resolutionRadii, detectorDistance, waveLength, 'DirectAxesP3P3rotation')
#logPredictedPattern(predictedPattern_P3P3, 'P3P3rotation')
