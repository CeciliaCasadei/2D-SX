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
import time
import joblib

def myParallelFunction(myObject, runNumber, detectorDistance, pixelSize, 
                       radialTolerance, pixelTolerance, azimuthTolerance, 
                       minNofPeaksPerLattice, maxNofPeaksPerImage, 
                       referenceCellSize):
    if myObject.selectionFlag == 1 and myObject.runNumber == runNumber and myObject.nPeaks > 0:
        myObject.indexingFunction(detectorDistance, pixelSize, 
                                  radialTolerance, pixelTolerance, azimuthTolerance,
                                  minNofPeaksPerLattice, maxNofPeaksPerImage, 
                                  referenceCellSize)

def latticeIndexingFunction(myArguments):
    
    resolutionRadii = [50, 10, 7] 
    
    # DEFAULTS
    minNofPeaksPerLattice = 16
    maxNofPeaksPerImage = 250
    
    # READ INPUTS
    string1 = 'Usage: python latticeIndexing.py --referenceCellSize <referenceCellSize> --runNumber <runNumber> --detectorDistance <detectorDistance>'
    string2 = ' --pixelSize <pixelSize> --radialTolerance <radialTolerance> --pixelTolerance <pixelTolerance> --azimuthTolerance <azimuthTolerance>' 
    string3 = ' --minNofPeaksPerLattice <minNofPeaksPerLattice> --maxNofPeaksPerImage <maxNofPeaksPerImage> '  
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["referenceCellSize=",
                                                                 "runNumber=", "detectorDistance=", "pixelSize=", 
                                                                 "radialTolerance=", "pixelTolerance=", "azimuthTolerance=", 
                                                                 "minNofPeaksPerLattice=", "maxNofPeaksPerImage="])
    except getopt.GetoptError:
        print string1 + string2 + string3
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print string1 + string2 + string3
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
        
    
    
    if os.path.exists('./Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber)):
        
        # LOAD diffractionImage OBJECTS
        fRead = open('./Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber), 'rb')
        imageObjects = pickle.load(fRead) 
        fRead.close()
        
        # LATTICE INDEXING            
        startTime = time.time()
        print 'Indexing started.'
        
        joblib.Parallel(n_jobs=1)(joblib.delayed(myParallelFunction)
                       (imageObjects['%s'%i], runNumber, detectorDistance, pixelSize, 
                        radialTolerance, pixelTolerance, azimuthTolerance, 
                        minNofPeaksPerLattice, maxNofPeaksPerImage, 
                        referenceCellSize) 
                        for i,j in imageObjects.items())
        
        # MERGE LATTICE DICTIONARIES TO ONE DICTIONARY
        allLatticesDictionary = {}
        nProcessedImages = 0
        for i, j in imageObjects.items():
            if imageObjects['%s' %i].selectionFlag == 1 and imageObjects['%s' %i].runNumber == runNumber: 
                nProcessedImages = nProcessedImages + 1
                myPath = './Output_r%s/LatticeIndexing/latticeDictionary_r%s_image_%s.pkl'%(runNumber, imageObjects['%s' %i].runNumber, imageObjects['%s' %i].imageNumber)
                if os.path.exists(myPath):  
                    openPkl = open(myPath, 'rb')
                    singleImageLatticesDictionary = pickle.load(openPkl)
                    openPkl.close()
                    allLatticesDictionary = dict(allLatticesDictionary, **singleImageLatticesDictionary)        
        if nProcessedImages == 0:
            print 'WARNING: No images were processed.'                 
        allLatticesDictionaryFile = open('./Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.pkl'%(runNumber, runNumber), 'wb')
        pickle.dump(allLatticesDictionary, allLatticesDictionaryFile)
        allLatticesDictionaryFile.close()

        # REPORT RUNTIME
        runTime = time.time() - startTime
        print 'Indexing took %.1f s'%runTime            

        # REPORT INDEXING RESULTS
        print '\nSTORING INDEXING LOG FOR RUN %s IN FILE ./Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.txt\n'%(runNumber, runNumber, runNumber)
        logFile = './Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.txt'%(runNumber, runNumber)  
        f = open(logFile, 'w') 
        f.write('INDEXING RESULTS FOR RUN %s\n'%runNumber)
        k = 0
        for i, j in allLatticesDictionary.items():
            k = k + 1
        f.write('Number of dictionary items (lattices): %d\n'%k)
        for i, j in sorted(allLatticesDictionary.items()):
            f.write('\nFile name: %s\n'%j.fileName)
            f.write('Run number: %s\n'%j.runNumber)
            f.write('Tilt angle: %s degrees\n'%j.tiltAngle)                
            f.write('Image number: %s\n'%j.imageNumber)
            f.write('Lattice number (in image): %s\n'%j.latticeNumberInImage)                                         
            f.write('Wavelength: %.3f A\n'%j.wavelength)
            f.write('Rotation angle: %.3f radians\n'%j.inPlaneRotation)
            f.write('Pixel size: %f m\n'%j.pixelSize)
            f.write('Detector distance: %.3f m\n'%j.detectorDistance)
            f.write('Number of matched peaks: %d\n'%j.nMatchedPeaks)
            f.write('Reference cell size: %.3f A\n'%j.referenceCellSize)
            f.write('Indexed peaks table:\n')
            f.write('    h    k    observed radius    observed azimuth    observed I    radial delta    azimuth delta    obs peak n    predicted peak n    predicted radius    predicted azimuth\n')
            for z in range(0, j.nMatchedPeaks):
                f.write('%5d%5d%19.3f%20.3f%14.2f%16.3f%17.3f%14d%20d%20.3f%21.3f\n'
                      %(j.indexedPeaksTable[z,0], j.indexedPeaksTable[z,1],
                        j.indexedPeaksTable[z,2], j.indexedPeaksTable[z,3], j.indexedPeaksTable[z,4],
                        j.indexedPeaksTable[z,5], j.indexedPeaksTable[z,6], j.indexedPeaksTable[z,7],
                        j.indexedPeaksTable[z,8], j.indexedPeaksTable[z,9], j.indexedPeaksTable[z,10]))
        f.close()
        
        # PLOT INDEXING RESULTS
        startPlotTime = time.time()
        print 'Plotting indexed experimental peaks.'
        for i, j in sorted(imageObjects.items()):
            if imageObjects['%s' %i].selectionFlag == 1 and imageObjects['%s' %i].runNumber == runNumber and imageObjects['%s' %i].nPeaks > 0: 
                imageObjects['%s' %i].plotIndexedExperimentalPeaks(detectorDistance, pixelSize, resolutionRadii)
        plotRunTime = time.time() - startPlotTime
        print 'Plotting took %.1f s'%plotRunTime 
        
        # REMOVE INDIVIDUAL LATTICE DICTIONARIES
        os.system('rm ./Output_r%s/LatticeIndexing/latticeDictionary*.pkl'%runNumber)
  
if __name__ == "__main__":
    print "\n**** CALLING latticeIndexing ****"
    latticeIndexingFunction(sys.argv[1:])