# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot
import math
import h5py
import unassembledMatching

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
        
        index = index + 1
    
    return predictedPattern    

### EXTRACT GEOMETRY ###
geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5' 
geometryData = h5py.File(geometryFile, 'r')
xGeometry = geometryData['/x']   ### float32 ###
xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
yGeometry = geometryData['/y']   ### float32 ###
yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)

# START
cellSize = 62.45
resolutionLimit = 7
hmax = 100
kmax = 100

directVectors = numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]])
identityMatrix = numpy.matrix([[1, 0],[0, 1]])
reciprocalLattice_I = buildReciprocalLatticeFunction(cellSize, identityMatrix, directVectors, hmax, kmax, resolutionLimit)

waveLength = 1.486567
waveVector = 2*numpy.pi / waveLength
tiltAngle = 15.0 * numpy.pi / 180
detectorDistance = 0.285
pixelSize = 0.000110

rodIndices = [[1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                                    
for indices in rodIndices:
    print indices
    hRod = indices[0]
    kRod = indices[1]

    qRods = []
    qRods_all = []
    
    for i in numpy.linspace(0.0, 360.0, num=1440):
        print i
        inPlaneRotation = i * numpy.pi / 180  
        
        predictedPattern_I = calculatePredictedPatternFunction(reciprocalLattice_I, 
                                                               inPlaneRotation, 
                                                               waveVector, 
                                                               tiltAngle, 
                                                               detectorDistance, 
                                                               pixelSize)
        
        for spot in predictedPattern_I:
            h = spot[0]
            k = spot[1]
            if (h == hRod and k == kRod) or (h == -hRod-kRod and k == hRod) or (h == kRod and k == -hRod-kRod):
                
                detectorX = pixelSize * spot[10] * math.cos(spot[8]) ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                detectorY = pixelSize * spot[10] * math.sin(spot[8]) ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                
                ### MATCH PREDICTED SPOT POSITION TO POSITION IN UNASSEMBLED MATRIX VIA GEOMETRY ###
                xIndices = numpy.argwhere(abs(xGeometry - detectorX) <= 0.000055)
                xIndices = numpy.asarray(xIndices, dtype = numpy.int16)
                yIndices = numpy.argwhere(abs(yGeometry - detectorY) <= 0.000055)
                yIndices = numpy.asarray(yIndices, dtype = numpy.int16)
                successFlag, i_good, j_good = unassembledMatching.unassembledMatching(xIndices, yIndices) 
                if successFlag == 1:
                    qRods.append(spot[11])
                qRods_all.append(spot[11])
                
            if (h == -hRod and k == -kRod) or (h == hRod+kRod and k == -hRod) or (h == -kRod and k == hRod+kRod):
                
                detectorX = pixelSize * spot[10] * math.cos(spot[8]) ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                detectorY = pixelSize * spot[10] * math.sin(spot[8]) ### POSITION IN m RELATIVE TO DIRECT BEAM ###
                
                ### MATCH PREDICTED SPOT POSITION TO POSITION IN UNASSEMBLED MATRIX VIA GEOMETRY ###
                xIndices = numpy.argwhere(abs(xGeometry - detectorX) <= 0.000055)
                xIndices = numpy.asarray(xIndices, dtype = numpy.int16)
                yIndices = numpy.argwhere(abs(yGeometry - detectorY) <= 0.000055)
                yIndices = numpy.asarray(yIndices, dtype = numpy.int16)
                successFlag, i_good, j_good = unassembledMatching.unassembledMatching(xIndices, yIndices)            
                if successFlag == 1:
                    qRods.append(-spot[11])
                qRods_all.append(-spot[11])
                    
    matplotlib.pyplot.figure(facecolor = 'w')
    (n, bins, patches) = matplotlib.pyplot.hist(qRods, bins=60)
    matplotlib.pyplot.savefig('/home/scratch/casadei_c/HISTOGRAM_Tilt15_rod_%d_%d.png'%(hRod, kRod))
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure(facecolor = 'w')
    (n, bins, patches) = matplotlib.pyplot.hist(qRods_all, bins=60)
    matplotlib.pyplot.savefig('/home/scratch/casadei_c/HISTOGRAM_NoGeometry_Tilt15_rod_%d_%d.png'%(hRod, kRod))
    matplotlib.pyplot.close()