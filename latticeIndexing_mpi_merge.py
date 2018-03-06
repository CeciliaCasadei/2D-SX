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
    outputPath = './Output_r%s/LatticeIndexing'%runNumber
    
    if os.path.exists('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber)):
        
        # LOAD diffractionImage OBJECTS
        fRead = open('%s/r%s_imagesDictionary.pkl'%(inputPath, runNumber), 'rb')
        imageObjects = pickle.load(fRead) 
        fRead.close()
        
        # MERGE LATTICE DICTIONARIES TO ONE DICTIONARY
        allLatticesDictionary = {}
        nProcessedImages = 0
        for i, j in imageObjects.items():
            if (imageObjects['%s' %i].selectionFlag == 1 and 
                imageObjects['%s' %i].runNumber == runNumber): 
                nProcessedImages = nProcessedImages + 1
                myPath = '%s/latticeDictionary_r%s_image_%s.pkl'%(outputPath, 
                                                                  imageObjects['%s' %i].runNumber, 
                                                                  imageObjects['%s' %i].imageNumber)
                if os.path.exists(myPath):  
                    openPkl = open(myPath, 'rb')
                    singleImageLatticesDictionary = pickle.load(openPkl)
                    openPkl.close()
                    allLatticesDictionary = dict(allLatticesDictionary, 
                                                 **singleImageLatticesDictionary)        
        if nProcessedImages == 0:
            print 'WARNING: No images were processed.'                 
        allLatticesDictionaryFile = open('%s/r%s_allLatticesDictionary.pkl'%(outputPath, 
                                                                             runNumber), 
                                                                             'wb')
        pickle.dump(allLatticesDictionary, allLatticesDictionaryFile)
        allLatticesDictionaryFile.close()

                

        # REPORT INDEXING RESULTS
        print '\nSTORING INDEXING LOG FOR RUN %s IN FILE ./Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.txt\n'%(runNumber, runNumber, runNumber)
        logFile = './Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.txt'%(runNumber, 
                                                                                runNumber)  
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
                      %(j.indexedPeaksTable[z,0], 
                        j.indexedPeaksTable[z,1],
                        j.indexedPeaksTable[z,2], 
                        j.indexedPeaksTable[z,3], 
                        j.indexedPeaksTable[z,4],
                        j.indexedPeaksTable[z,5], 
                        j.indexedPeaksTable[z,6], 
                        j.indexedPeaksTable[z,7],
                        j.indexedPeaksTable[z,8], 
                        j.indexedPeaksTable[z,9], 
                        j.indexedPeaksTable[z,10]))
        f.close()
        
        # PLOT INDEXING RESULTS
        plotFlag = 0
        if plotFlag == 1:
            print 'Plotting indexed experimental peaks.'
            for i, j in sorted(imageObjects.items()):
                if (imageObjects['%s' %i].selectionFlag == 1 and 
                    imageObjects['%s' %i].runNumber == runNumber and 
                    imageObjects['%s' %i].nPeaks > 0): 
                    imageObjects['%s' %i].plotIndexedExperimentalPeaks(detectorDistance, 
                                                                       pixelSize, 
                                                                       resolutionRadii)
            
        
        # REMOVE INDIVIDUAL LATTICE DICTIONARIES
        os.system('rm ./Output_r%s/LatticeIndexing/latticeDictionary*.pkl'%runNumber)
  
if __name__ == "__main__":
    mergeFunction(sys.argv[1:])