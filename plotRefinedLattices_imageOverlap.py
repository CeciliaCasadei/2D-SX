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
import time
import joblib

imageToPlot = '892'
def myParallelFunction(myObject, runNumber, resolutionRadii):
    if myObject.selectionFlag == 1 and myObject.runNumber == runNumber and myObject.nPeaks > 0:
        if myObject.imageNumber == imageToPlot:
            myObject.plotRefinedLattices_imageOverlap(resolutionRadii)

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
    
    imagesDictionaryFile = './Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber)
    if not os.path.exists(imagesDictionaryFile):
        print 'File %s not found.'%imagesDictionaryFile
    else:
        fRead = open(imagesDictionaryFile, 'rb')
        myData = pickle.load(fRead)                                            # Dictionary items are diffractionImage objects
        fRead.close()
    
        startTime = time.time()
        print 'Plotting started.'
        joblib.Parallel(n_jobs=1)(joblib.delayed(myParallelFunction)(myData['%s'%i], runNumber, resolutionRadii) for i,j in myData.items())
                            
        runTime = time.time() - startTime
        print 'Plotting took %.1f s.'%runTime

if __name__ == "__main__":
    print "\n**** CALLING plotRefinedLattices ****"
    plotRefinedLatticesFunction(sys.argv[1:]) 