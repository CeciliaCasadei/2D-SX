# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 15:10:58 2016
@author: casadei_c
LAUNCH PARALLEL INDEXING (IMAGES ASSIGNED TO 4 PROCESSORS)
AND PLOT INDEXING RESULTS.
"""

import sys
import getopt
import os
import pickle
from mpi4py import MPI



def mergeFunction(myArguments):
    
    resolutionRadii = [50, 10, 7] 
     
    # READ INPUTS
    string1 = 'Usage: python latticeIndexing.py'
    string2 = ' --runNumber <runNumber>'
    string3 = ' --detectorDistance <detectorDistance>'
    string4 = ' --pixelSize <pixelSize>' 
    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", 
                                                                 "detectorDistance=", 
                                                                 "pixelSize="])
    except getopt.GetoptError:
        print string1 + string2 + string3 + string4 
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print string1 + string2 + string3 + string4 
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--detectorDistance":
            detectorDistance = float(value)
        elif option == "--pixelSize":
            pixelSize = float(value)

        
    
    inputPath = './Output_r%s/ExtractExperimentalInfo'%runNumber

    
    if os.path.exists('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber)):
        
        # LOAD diffractionImage OBJECTS
        fRead = open('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber), 'rb')
        imageObjects = pickle.load(fRead) 
        fRead.close()
        
        # PLOT INDEXING RESULTS
        plotFlag = 1
        if plotFlag == 1:
            
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            numberOfCores = comm.Get_size()
            if rank == 0:
                print "I see %d cores are available"%(numberOfCores)            
            
            iImage = 0
            for i, j in sorted(imageObjects.items()):
                
                iImage += 1
                if iImage%numberOfCores != rank:
                    continue
            
                if (j.selectionFlag == 1 and 
                    j.runNumber == runNumber and 
                    j.nPeaks > 0): 
                    j.plotIndexedExperimentalPeaks(detectorDistance, 
                                                   pixelSize, 
                                                   resolutionRadii)
            
        
        # REMOVE INDIVIDUAL LATTICE DICTIONARIES
        #os.system('rm ./Output_r%s/LatticeIndexing/latticeDictionary*.pkl'%runNumber)
  
if __name__ == "__main__":
    mergeFunction(sys.argv[1:])