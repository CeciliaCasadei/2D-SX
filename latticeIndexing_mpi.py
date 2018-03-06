# -*- coding: utf-8 -*-
"""
Created on 05/03/2018
@author: casadei_c
LAUNCH PARALLEL INDEXING 
"""

import sys
import getopt
import os
import pickle
from mpi4py import MPI


def latticeIndexingFunction(myArguments):
    
    # DEFAULTS
    minNofPeaksPerLattice = 16
    maxNofPeaksPerImage = 250
    
    # READ INPUTS
    string1 = 'Usage: python latticeIndexing.py --referenceCellSize <referenceCellSize>'
    string2 = ' --runNumber <runNumber> --detectorDistance <detectorDistance>'
    string3 = ' --pixelSize <pixelSize> --radialTolerance <radialTolerance>' 
    string4 = ' --pixelTolerance <pixelTolerance> --azimuthTolerance <azimuthTolerance>'
    string5 = ' --minNofPeaksPerLattice <minNofPeaksPerLattice>'  
    string6 = ' --maxNofPeaksPerImage <maxNofPeaksPerImage>'
    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["referenceCellSize=",
                                                                 "runNumber=", 
                                                                 "detectorDistance=", 
                                                                 "pixelSize=", 
                                                                 "radialTolerance=", 
                                                                 "pixelTolerance=", 
                                                                 "azimuthTolerance=", 
                                                                 "minNofPeaksPerLattice=", 
                                                                 "maxNofPeaksPerImage="])
    except getopt.GetoptError:
        print string1 + string2 + string3 + string4 + string5 + string6
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print string1 + string2 + string3 + string4 + string5 + string6
            sys.exit()
        elif option == "--referenceCellSize":
            referenceCellSize = float(value)
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--detectorDistance":
            detectorDistance = float(value)
        elif option == "--pixelSize":
            pixelSize = float(value)
        elif option == "--radialTolerance":
            radialTolerance = float(value)    
        elif option == "--pixelTolerance":
            pixelTolerance = float(value)   
        elif option == "--azimuthTolerance":
            azimuthTolerance = float(value)  
        elif option == "--minNofPeaksPerLattice":
            minNofPeaksPerLattice = int(value)
        elif option == "--maxNofPeaksPerImage":
            maxNofPeaksPerImage = int(value)
        
    
    inputPath = './Output_r%s/ExtractExperimentalInfo'%runNumber
    if os.path.exists('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber)):
        
        # LOAD diffractionImage OBJECTS
        fRead = open('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber), 'rb')
        imageObjects = pickle.load(fRead) 
        fRead.close()
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        numberOfCores = comm.Get_size()
        if rank == 0:
            print "I see %d cores are available"%(numberOfCores)
        
        # Loop on all images
        iImage = 0
        for i,j in imageObjects.items():
            iImage += 1
            if iImage%numberOfCores != rank:
                continue

            if (imageObjects['%s'%i].selectionFlag == 1 and 
                imageObjects['%s'%i].runNumber == runNumber and 
                imageObjects['%s'%i].nPeaks > 0):
                    
                imageObjects['%s'%i].indexingFunction(detectorDistance, 
                                                      pixelSize, 
                                                      radialTolerance, 
                                                      pixelTolerance, 
                                                      azimuthTolerance,
                                                      minNofPeaksPerLattice, 
                                                      maxNofPeaksPerImage, 
                                                      referenceCellSize)

  
if __name__ == "__main__":
    latticeIndexingFunction(sys.argv[1:])