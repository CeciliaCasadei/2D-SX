# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 09:50:38 2016
@author: casadei_c
CALL diffractionImage CLASS FUNCTION plotRefinedLattices
"""
import sys
import getopt
import os
import pickle
from mpi4py import MPI

def parallelFunction(I, runNumber, resolutionRadii):
    if I.selectionFlag == 1 and I.runNumber == runNumber and I.nPeaks > 0:
        I.plotRefinedLattices(resolutionRadii)

def plotRefinedLatticesFunction(myArguments):
    resolutionRadii = [50, 10, 7]   
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python plotRefinedLattices.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRefinedLattices.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
    
    inputPath = './Output_r%s/ExtractExperimentalInfo'%runNumber
    imagesDictionaryFile = '%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber)
    
    if not os.path.exists(imagesDictionaryFile):
        print 'File %s not found.'%imagesDictionaryFile
    else:
        fRead = open(imagesDictionaryFile, 'rb')
        myData = pickle.load(fRead)                                            # Dictionary items are diffractionImage objects
        fRead.close()
        
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        numberOfCores = comm.Get_size()
        if rank == 0:
            print "I see %d cores are available"%(numberOfCores)
    
    
        iImg = 0
        for i,j in myData.items():
            iImg += 1
            if iImg%numberOfCores != rank:
                continue
        
            parallelFunction(j, runNumber, resolutionRadii)
                            
        

if __name__ == "__main__":
    #print "\n**** CALLING plotRefinedLattices ****"
    plotRefinedLatticesFunction(sys.argv[1:]) 