# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os

def scaling_applyScalesFunction(myArguments):
    
    # DEFAULTS
    outputFolder = ''
    
    # READ INPUTS 
    inputString = '--runNumber <runNumber> --outputFolder <outputFolder>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                             ["runNumber=", 
                                              "outputFolder="])
    except getopt.GetoptError:
        print 'Usage: python scaling_applyScales.py %s'%inputString
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_applyScales.py %s'%inputString
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--outputFolder":
            outputFolder = value

    if outputFolder == '':
        outputFolder = './Output_r%s/transformAndScale'%runNumber
    
    matricesListFolder = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, 
                                                                 runNumber)
    spotsMatricesList = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                                     %(matricesListFolder, 
                                       runNumber)) # h_transformed 
                                                   # k_transformed 
                                                   # qRod 
                                                   # I 
                                                   # flag 
                                                   # i_unassembled 
                                                   # j_unassembled
    latticeScales = joblib.load('%s/r%s-finalScales/r%s-finalScales-normalized.jbl'
                                 %(outputFolder, runNumber, runNumber))
    nLattices = len(latticeScales) 
    
    scaledSpotMatricesList = []
    for lattice in range(0, nLattices):
        latticeScale = latticeScales[lattice]
        spotsMatrix = spotsMatricesList[lattice]    
        print 'Lattice: %d - Scale: %.3f'%(lattice, latticeScale)
        if numpy.isnan(latticeScale):
            flag = 0
            scaledSpotsMatrix = []
            for spot in spotsMatrix:
                scaledSpot = [spot[0], 
                              spot[1], 
                              spot[2], 
                              spot[3], 
                              flag, 
                              spot[5], 
                              spot[6], 
                              latticeScale]
                scaledSpotsMatrix.append(scaledSpot)
        else:
            flag = 1
            scaledSpotsMatrix = []
            for spot in spotsMatrix:
                scaledI = spot[3] * latticeScale
                scaledSpot = [spot[0],       # h_transformed
                              spot[1],       # k_transformed
                              spot[2],       # qRod
                              scaledI,       # Iscaled
                              flag,          # flag
                              spot[5],       # i_unassembled
                              spot[6],       # j_unassembled
                              latticeScale]  # scale
                scaledSpotsMatrix.append(scaledSpot)
        scaledSpotsMatrix = numpy.asarray(scaledSpotsMatrix)                         
        scaledSpotMatricesList.append(scaledSpotsMatrix)
           
    outputPath = '%s/spotsMatricesList-Scaled-r%s'%(outputFolder, runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(scaledSpotMatricesList, 
                '%s/r%s_scaledSpotsMatricesList.jbl'%(outputPath, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING scaling_applyScales ****"
    scaling_applyScalesFunction(sys.argv[1:])    