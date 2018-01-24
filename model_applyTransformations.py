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


def model_applyTransformationsFunction(myArguments):
    
    inputStr = '--runNumber <runNumber>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python model_applyTransformations.py %s'%inputStr
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyTransformations.py %s'%inputStr
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    inFolder_path = './Output_r%s/transformAndScale'%runNumber
    inputFolder = '%s/spotsMatricesList-r%s'%(inFolder_path, runNumber)
    outFolder_path = './Output_runMergingVsModel'
    outputFolder = '%s/transformAndScaleToModel_r%s'%(outFolder_path,
                                                      runNumber)

    spotsMatricesList = joblib.load('%s/r%s_spotsMatricesList.jbl'%(inputFolder, 
                                                                    runNumber))# h 
                                                                               # k 
                                                                               # qRod 
                                                                               # I 
                                                                               # flag=1 
                                                                               # i_unassembled 
                                                                               # j_unassembled
    latticeOrientations = joblib.load('%s/r%s_orientationsVsModel.jbl'%(outputFolder, 
                                                                        runNumber))
    nLattices = len(latticeOrientations)
    
    transformedSpotMatricesList = []
    for lattice in range(0, nLattices):
        latticeOrientation = latticeOrientations[lattice]
        spotsMatrix = spotsMatricesList[lattice]                                                                              
        print 'Lattice: %d - Orientation Vs model: %s'%(lattice, 
                                                        latticeOrientation)
        if numpy.isnan(latticeOrientation):
            flag = 0
            transformedSpotsMatrix = []
            for spot in spotsMatrix:
                transformedSpot = [spot[0], 
                                   spot[1], 
                                   spot[2], 
                                   spot[3], 
                                   flag, 
                                   spot[5], 
                                   spot[6]]
                transformedSpotsMatrix.append(transformedSpot)
        else:
            flag = 1
            transformedSpotsMatrix = []
            latticeOrientationMatrix = indexToMatrix(latticeOrientation)
            for spot in spotsMatrix:
                h = spot[0]
                k = spot[1]
                indices = numpy.matrix('%d; %d'%(h, k)) #(2x1)
                transformedIndices = latticeOrientationMatrix*indices
                h_t = transformedIndices[0, 0]
                k_t = transformedIndices[1, 0]
                transformedSpot = [h_t, 
                                   k_t, 
                                   spot[2], 
                                   spot[3], 
                                   flag, 
                                   spot[5], 
                                   spot[6]]
                transformedSpotsMatrix.append(transformedSpot)
        transformedSpotsMatrix = numpy.asarray(transformedSpotsMatrix, 
                                               dtype = numpy.float32)                                
        transformedSpotMatricesList.append(transformedSpotsMatrix)
    
    outputPath = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, 
                                                         runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(transformedSpotMatricesList, 
                '%s/r%s_transformedSpotsMatricesList.jbl'%(outputPath, 
                                                           runNumber))

if __name__ == "__main__":
    print "\n**** CALLING model_applyTransformations ****"
    model_applyTransformationsFunction(sys.argv[1:])    