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
        print 'Usage: python model_applyTransformations.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyTransformations.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    inputFolder = './Output_r%s/transformAndScale'%runNumber
    outputFolder = './Output_runMerging/transformAndScaleToModel_r%s'%runNumber


    spotsMatricesList = joblib.load('%s/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'%(inputFolder, runNumber, runNumber))
    latticeOrientations = joblib.load('%s/r%s_orientationsVsModel.jbl'%(outputFolder, runNumber))
    nLattices = len(latticeOrientations)  # equal to len(spotsMatricesList)
    
    transformedSpotMatricesList = []
    for lattice in range(0, nLattices):
        latticeOrientation = latticeOrientations[lattice]
        spotsMatrix = spotsMatricesList[lattice] # n h k qRod I Icorrected
        print 'Lattice: %d - Orientation: %s - spotsMatrix: (%d, %d)'%(lattice, latticeOrientation, spotsMatrix.shape[0], spotsMatrix.shape[1])
        if numpy.isnan(latticeOrientation):
            transformedSpotsMatrix = spotsMatrix
        else:
            transformedSpotsMatrix = []
            latticeOrientationMatrix = indexToMatrix(latticeOrientation)
            for spot in spotsMatrix:
                h = spot[1]
                k = spot[2]
                indices = numpy.matrix('%d; %d'%(h, k)) #(2x1)
                transformedIndices = latticeOrientationMatrix*indices
                h_t = transformedIndices[0, 0]
                k_t = transformedIndices[1, 0]
                transformedSpot = [spot[0], spot[1], spot[2], spot[3], spot[4], spot[5], h_t, k_t]
                transformedSpotsMatrix.append(transformedSpot)
        transformedSpotsMatrix = numpy.asarray(transformedSpotsMatrix) # n h k qRod I Icorrected h_transformed k_transformed
        transformedSpotMatricesList.append(transformedSpotsMatrix)
        print 'Lattice: %d - Orientation: %s - transformedSpotsMatrix: (%d, %d)'%(lattice, latticeOrientation, transformedSpotsMatrix.shape[0], transformedSpotsMatrix.shape[1])
    
    outputPath = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(transformedSpotMatricesList, '%s/r%s_transformedSpotsMatricesList.jbl'%(outputPath, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING transform_applyTransformations ****"
    transform_applyTransformationsFunction(sys.argv[1:])    