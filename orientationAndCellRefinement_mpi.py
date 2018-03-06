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
import joblib
from mpi4py import MPI



def parallelFunction(L, 
                     sizeRefinementSteps, 
                     orientationRefinementSteps, 
                     folderName, 
                     resolutionLimit):
                         
    refinedParameters = L.refineCellSizeAndOrientation(sizeRefinementSteps, 
                                                       orientationRefinementSteps, 
                                                       folderName, 
                                                       resolutionLimit)    
                                                       
    L.setRefinedInPlaneOrientation(refinedParameters['refinedOrientation'])   # Set Lattice object attribute: refinedInPlaneOrientation
    L.setRefinedCellSize(refinedParameters['refinedCellSize'])                # Set Lattice object attribute: refinedCellSize
    L.setLatticeErrorMatrix(refinedParameters['latticeErrorMatrix'])          # Set Lattice object attribute: latticeErrorMatrix
    L.setRefinedPattern(refinedParameters['refinedPredictedPattern'])         # Set Lattice object attribute: refinedPredictedPattern
    L.setLatticeError(refinedParameters['averageError'])                      # Set Lattice object attribute: avgLatticeError
    joblib.dump(L, '%s/r%s_Image%s_Lattice%s.jbl'%(folderName, 
                                                   L.runNumber, 
                                                   L.imageNumber, 
                                                   L.latticeNumberInImage))
    L.logRefinedPatternInfo('%s/RefinedPredictedPatterns'%folderName)
    L.plotLatticeErrorMatrix('%s/Figures/r%s_Image%s_Lattice%s.png'%(folderName, 
                                                                     L.runNumber, 
                                                                     L.imageNumber, 
                                                                     L.latticeNumberInImage))
       
    
def orientationAndCellRefinementFunction(myArguments): 
        
    # DEFAULTS
    nSizeRefSteps = 21                  # int
    nOrientationRefSteps = 21           # int
    widthSizeRefSteps = 0.004           # fraction
    widthOrientationRefSteps = 0.2      # degrees
    resolutionLimit = 5                 # A
    
    # READ INPUTS
    str1 = 'Usage: python orientationAndCellRefinement.py'
    str2 = ' --referenceCellSize <referenceCellSize>'
    str3 = ' --runNumber <runNumber>'
    str4 = ' --nSizeRefSteps <nSizeRefSteps>'
    str5 = ' --nOrientationRefSteps <nOrientationRefSteps>'  
    str6 = ' --widthSizeRefSteps <widthSizeRefSteps>'
    str7 = ' --widthOrientationRefSteps <widthOrientationRefSteps>'
    str8 = ' --resolutionLimit <resolutionLimit>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["referenceCellSize=", 
                                                                 "runNumber=", 
                                                                 "nSizeRefSteps=", 
                                                                 "nOrientationRefSteps=", 
                                                                 "widthSizeRefSteps=", 
                                                                 "widthOrientationRefSteps=", 
                                                                 "resolutionLimit="])
    except getopt.GetoptError:
        print str1 + str2 + str3 + str4 + str5 + str6 + str7 + str8
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print str1 + str2 + str3 + str4 + str5 + str6 + str7 + str8
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
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)                                     # A
    
    latticeFolder = '/Output_r%s/LatticeIndexing'%runNumber
    refinementFolder = './Output_r%s/OrientationAndCellRefinement'%runNumber       
    
    # DEFINE SIZE SEARCH LIST 
    nSizeRefStepsOverTwo = nSizeRefSteps / 2
    sizeRefWidth = nSizeRefStepsOverTwo*widthSizeRefSteps
    sizeRefinementRescalings = numpy.linspace(1-sizeRefWidth, 
                                              1+sizeRefWidth, 
                                              nSizeRefSteps)
    sizeRefinementSteps = referenceCellSize * sizeRefinementRescalings    

    # DEFINE ORIENTATION SEARCH LIST
    nOrientationRefStepsOverTwo = nOrientationRefSteps / 2
    widthOrientationRefSteps = widthOrientationRefSteps/180*numpy.pi
    orientationRefWidth = nOrientationRefStepsOverTwo*widthOrientationRefSteps
    orientationRefinementSteps = numpy.linspace(-orientationRefWidth,
                                                +orientationRefWidth,
                                                nOrientationRefSteps) 
                                                
    # LOAD LATTICES
    latticeDictionaryFile = '.%s/r%s_allLatticesDictionary.pkl'%(latticeFolder, 
                                                                 runNumber)      
    if not os.path.exists(latticeDictionaryFile):
        print 'File %s not found.'%latticeDictionaryFile
        return          
    fRead = open(latticeDictionaryFile, 'rb')
    myData = pickle.load(fRead)                                                
    fRead.close()
    
    # ORIENTATION AND CELL SIZE REFINEMENT
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numberOfCores = comm.Get_size()
    if rank == 0:
        print "I see %d cores are available"%(numberOfCores)
        
    iLattice = 0
    for i,j in myData.items():
        iLattice += 1
        if iLattice%numberOfCores != rank:
            continue
        
        parallelFunction(j, 
                         sizeRefinementSteps, 
                         orientationRefinementSteps, 
                         refinementFolder, 
                         resolutionLimit)
        



if __name__ == "__main__":
    #print "\n**** CALLING orientationAndCellRefinement ****"
    orientationAndCellRefinementFunction(sys.argv[1:])