# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os

def scaling_applyScalesFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python scaling_applyScales.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_applyScales.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_r%s/transformAndScale'%runNumber
    spotsMatricesList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    latticeScales = joblib.load('%s/r%s-finalScales.jbl'%(outputFolder, runNumber))
    nLattices = len(latticeScales)  # equal to len(spotsMatricesList)
    
    scaledSpotMatricesList = []
    for lattice in range(0, nLattices):
        latticeScale = latticeScales[lattice]
        spotsMatrix = spotsMatricesList[lattice] # n h k qRod I Icorrected h_tranformed k_transformed
        print 'Lattice: %d - Scale: %.3f - spotsMatrix: (%d, %d)'%(lattice, latticeScale, spotsMatrix.shape[0], spotsMatrix.shape[1])
        if numpy.isnan(latticeScale):
            scaledSpotsMatrix = spotsMatrix
        else:
            scaledSpotsMatrix = []
            for spot in spotsMatrix:
                scaledI = spot[5] * latticeScale
                scaledSpot = [spot[0], spot[1], spot[2], spot[3], spot[4], spot[5], spot[6], spot[7], scaledI]
                scaledSpotsMatrix.append(scaledSpot)
        scaledSpotsMatrix = numpy.asarray(scaledSpotsMatrix) # n h k qRod I Icorrected h_transformed k_transformed
        scaledSpotMatricesList.append(scaledSpotsMatrix)
        print 'Lattice: %d - Scale: %.3f - scaledSpotMatricesList: (%d, %d)'%(lattice, latticeScale, scaledSpotsMatrix.shape[0], scaledSpotsMatrix.shape[1])
    
    outputPath = '%s/spotsMatricesList-Scaled-r%s'%(outputFolder, runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(scaledSpotMatricesList, '%s/r%s_scaledSpotsMatricesList.jbl'%(outputPath, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING scaling_applyScales ****"
    scaling_applyScalesFunction(sys.argv[1:])    