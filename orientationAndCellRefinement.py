# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 17:06:00 2016
@author: casadei_c
READ ARGUMENTS FOR ORIENTATION AND CELL SIZE REFINEMENT.
BUILD CELL SIZE AND IN-PLANE ORIENTATION SEARCH LISTS.
FOR EACH TRIAL CELL SIZE, CALCULATE THE ASSOCIATED RECIPROCAL LATTICE.
FOR EACH LATTICE, REFINE CELL SIZE AND IN PLANE ORIENTATION
BY CALLING Lattice CLASS FUNCTION refineCellSizeAndOrientation
(PARALLELIZED ON LATTICES).
"""
import sys
import getopt
import os
import numpy
import pickle
import time
import joblib

import buildReciprocalLattice

def myParallelFunction(myObject, sizeRefinementSteps, orientationRefinementSteps, folderName, RLfolderName):
    refinedParameters = myObject.refineCellSizeAndOrientation(sizeRefinementSteps, 
                                                              orientationRefinementSteps, 
                                                              folderName, RLfolderName)    
    myObject.setRefinedInPlaneOrientation(refinedParameters['refinedOrientation'])   # Set Lattice object attribute: refinedInPlaneOrientation
    myObject.setRefinedCellSize(refinedParameters['refinedCellSize'])                # Set Lattice object attribute: refinedCellSize
    myObject.setLatticeErrorMatrix(refinedParameters['latticeErrorMatrix'])          # Set Lattice object attribute: latticeErrorMatrix
    myObject.setRefinedPattern(refinedParameters['refinedPredictedPattern'])         # Set Lattice object attribute: refinedPredictedPattern
    myObject.setLatticeError(refinedParameters['averageError'])                      # Set Lattice object attribute: avgLatticeError
    joblib.dump(myObject, '%s/r%s_Image%s_Lattice%s.jbl'
                           %(folderName, myObject.runNumber, myObject.imageNumber, myObject.latticeNumberInImage))
    myObject.logRefinedPatternInfo('%s/RefinedPredictedPatterns'%folderName)
    
def orientationAndCellRefinementFunction(myArguments): 
        
    # DEFAULTS
    nSizeRefSteps = 21                  # int
    nOrientationRefSteps = 21           # int
    widthSizeRefSteps = 0.004           # fraction
    widthOrientationRefSteps = 0.2      # degrees
    resolutionLimit = 5                 # A
    
    # READ INPUTS
    string1 = 'Usage: python orientationAndCellRefinement.py --referenceCellSize <referenceCellSize> --runNumber <runNumber> --nSizeRefSteps <nSizeRefSteps>'
    string2 = ' --nOrientationRefSteps <nOrientationRefSteps> --widthSizeRefSteps <widthSizeRefSteps> --widthOrientationRefSteps <widthOrientationRefSteps>'    
    string3 = ' --hmax <hmax> --kmax <kmax> --resolutionLimit <resolutionLimit>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["referenceCellSize=", "runNumber=", 
                                                                 "nSizeRefSteps=", "nOrientationRefSteps=", 
                                                                 "widthSizeRefSteps=", "widthOrientationRefSteps=", 
                                                                 "hmax=", "kmax=",
                                                                 "resolutionLimit="])
    except getopt.GetoptError:
        print string1 + string2 + string3
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print string1 + string2 + string3
            sys.exit()
        elif option == "--referenceCellSize":
            referenceCellSize = float(value)                                   # A
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--nSizeRefSteps":
            nSizeRefSteps = int(value)
        elif option == "--nOrientationRefSteps":
            nOrientationRefSteps = int(value)
        elif option == "--widthSizeRefSteps":
            widthSizeRefSteps = float(value)    
        elif option == "--widthOrientationRefSteps":
            widthOrientationRefSteps = float(value)   
        elif option == "--hmax":
            hmax = int(value)  
        elif option == "--kmax":
            kmax = int(value)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)                                     # A
            
    # FOLDERS
    refinementFolder = './Output_r%s/OrientationAndCellRefinement'%runNumber
    if not os.path.exists(refinementFolder):
        os.mkdir(refinementFolder)
    refinedPatternsFolder = '%s/RefinedPredictedPatterns'%refinementFolder
    if not os.path.exists(refinedPatternsFolder):
        os.mkdir(refinedPatternsFolder)
    reciprocalLatticesFolder = '%s/ReciprocalLattices'%refinementFolder
    if not os.path.exists(reciprocalLatticesFolder):
        os.mkdir(reciprocalLatticesFolder)
    figuresFolder = '%s/Figures'%refinementFolder
    if not os.path.exists(figuresFolder):
        os.mkdir(figuresFolder)
    
    
    # DEFINE SIZE SEARCH LIST AND BUILD RECIPROCAL LATTICES
    nSizeRefStepsOverTwo = nSizeRefSteps / 2
    sizeRefinementRescalings = numpy.linspace(1-nSizeRefStepsOverTwo*widthSizeRefSteps, 
                                              1+nSizeRefStepsOverTwo*widthSizeRefSteps, 
                                              nSizeRefSteps)
    sizeRefinementSteps = referenceCellSize * sizeRefinementRescalings    
    for sizeRefinementStep in sizeRefinementSteps:
        reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(sizeRefinementStep, hmax, kmax, resolutionLimit)
        f = open('./Output_r%s/OrientationAndCellRefinement/ReciprocalLattices/reciprocalLattice_cellSize_%.3f.pkl'%(runNumber, sizeRefinementStep), 'wb') 
        pickle.dump(reciprocalLattice, f)
        f.close()
    
    # DEFINE ORIENTATION SEARCH LIST
    nOrientationRefStepsOverTwo = nOrientationRefSteps / 2
    widthOrientationRefSteps = widthOrientationRefSteps/180*numpy.pi
    orientationRefinementSteps = numpy.linspace(-nOrientationRefStepsOverTwo*widthOrientationRefSteps,
                                                +nOrientationRefStepsOverTwo*widthOrientationRefSteps,
                                                nOrientationRefSteps) 
                                                
    # LOAD LATTICES
    latticeDictionaryFile = './Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.pkl'%(runNumber, runNumber)      
    if not os.path.exists(latticeDictionaryFile):
        print 'File %s not found.'%latticeDictionaryFile
        return          
    fRead = open(latticeDictionaryFile, 'rb')
    myData = pickle.load(fRead)                                                ### Dictionary items are Lattice objects ###
    fRead.close()
    
    # ORIENTATION AND CELL SIZE REFINEMENT
    startTime = time.time()
    print 'Refinement started.'
    joblib.Parallel(n_jobs=-1)(joblib.delayed(myParallelFunction)
                                             (myData['%s'%i], sizeRefinementSteps, orientationRefinementSteps, 
                                             refinementFolder, reciprocalLatticesFolder) 
                                             for i,j in myData.items())
    
    # MERGE REFINED LATTICES TO ONE DICTIONARY
    allRefinedLatticesDictionary = {}
    for i, j in myData.items():                                                                  
        myPath = '%s/r%s_Image%s_Lattice%s.jbl'%(refinementFolder, 
                                                 myData['%s' %i].runNumber, 
                                                 myData['%s' %i].imageNumber, 
                                                 myData['%s' %i].latticeNumberInImage)
        if os.path.exists(myPath):                 
            myRefinedLattice = joblib.load(myPath)
            allRefinedLatticesDictionary['%s_Lattice_%d'%(myRefinedLattice.fileName, 
                                                          myRefinedLattice.latticeNumberInImage)] = myRefinedLattice                  
        else:
            print 'File %s not found.'%myPath                
    allLatticesDictionaryFile = open('%s/r%s_refinedLatticesDictionary.pkl'%(refinementFolder, runNumber), 'wb')
    pickle.dump(allRefinedLatticesDictionary, allLatticesDictionaryFile)
    allLatticesDictionaryFile.close()
    
    # LOG REFINEMENT RESULTS    
    fSummary = open('%s/r%s_Summary.txt'%(refinementFolder, runNumber), 'w')
    for i,j in sorted(allRefinedLatticesDictionary.items()):
        fSummary.write('Run: %s'%j.runNumber)
        fSummary.write('\nImage n: %s'%j.imageNumber)
        fSummary.write('\nLattice n: %s'%j.latticeNumberInImage)
        fSummary.write('\nTilt angle: %s degrees'%j.tiltAngle)
        fSummary.write('\nOriginal cell size: %.2f A'%j.referenceCellSize)
        fSummary.write('\nRefined cell size: %.2f A'%j.refinedCellSize)
        fSummary.write('\nOriginal in-plane orientation: %.3f radians'%j.inPlaneRotation)
        fSummary.write('\nRefined in-plane orientation: %.3f radians'%j.refinedInPlaneOrientation)
        fSummary.write('\n\n')
    fSummary.close()
    runTime = time.time() - startTime
    print 'Cell size and orientation refinement took %.2f s'%runTime 
    
    # REMOVE INDIVIDUAL LATTICE FILES
    os.system('rm %s/*.jbl'%refinementFolder)
    os.system('rm %s/*.npy'%refinementFolder)

    # PLOT ERROR MATRICES
    print 'Plotting error matrices'     
    for i,j in allRefinedLatticesDictionary.items():
        j.plotLatticeErrorMatrix('%s/r%s_Image%s_Lattice%s.png'
                                  %(figuresFolder, j.runNumber, j.imageNumber, j.latticeNumberInImage))
        

if __name__ == "__main__":
    print "\n**** CALLING orientationAndCellRefinement ****"
    orientationAndCellRefinementFunction(sys.argv[1:])