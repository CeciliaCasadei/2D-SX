# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os


def indexToMatrix(index):
    identity = numpy.matrix([[1, 0],[0, 1]])
    inversion = numpy.matrix([[-1, 0],[0, -1]])
    permutation = numpy.matrix([[0, 1],[1, 0]])
    inversion_permutation = numpy.matrix([[0, -1],[-1, 0]])    
    if index == 0:
        matrix = identity
    elif index == 1:
        matrix = inversion
    elif index == 2:
        matrix = permutation
    elif index == 3:
        matrix = inversion_permutation
    else:
        matrix = numpy.isnan
    return matrix


def transform_applyTransformationsFunction(myArguments):
    
    runNumber = ''
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python transform_applyTransformations.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_applyTransformations.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_r%s/transformAndScale'%runNumber


    spotsMatricesList = joblib.load('%s/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    latticeOrientations = joblib.load('%s/r%s-finalOrientations.jbl'%(outputFolder, runNumber))
    nLattices = len(latticeOrientations)  # equal to len(spotsMatricesList)
    
    transformedSpotMatricesList = []
    for lattice in range(0, nLattices):
        latticeOrientation = latticeOrientations[lattice]
        spotsMatrix = spotsMatricesList[lattice] # h k qRod I 
        print 'Lattice: %d - Orientation: %s'%(lattice, latticeOrientation)
        if numpy.isnan(latticeOrientation):
            flag = 0
            transformedSpotsMatrix = []
            for spot in spotsMatrix:
                h = spot[0]
                k = spot[1]
                q = spot[2]
                I = spot[3]
                transformedSpot = [h, k, q, I, flag]
                transformedSpotsMatrix.append(transformedSpot)
        else:
            flag = 1
            transformedSpotsMatrix = []
            latticeOrientationMatrix = indexToMatrix(latticeOrientation)
            for spot in spotsMatrix:
                h = spot[0]
                k = spot[1]
                q = spot[2]
                I = spot[3]
                indices = numpy.matrix('%d; %d'%(h, k)) #(2x1)
                transformedIndices = latticeOrientationMatrix*indices
                h_t = transformedIndices[0, 0]
                k_t = transformedIndices[1, 0]
                transformedSpot = [h_t, k_t, q, I, flag]
                transformedSpotsMatrix.append(transformedSpot)
                
        transformedSpotsMatrix = numpy.asarray(transformedSpotsMatrix) # h_transformed k_transformed qRod I flag
        transformedSpotMatricesList.append(transformedSpotsMatrix)
        
    
    outputPath = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(transformedSpotMatricesList, '%s/r%s_transformedSpotsMatricesList.jbl'%(outputPath, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING transform_applyTransformations ****"
    transform_applyTransformationsFunction(sys.argv[1:])    