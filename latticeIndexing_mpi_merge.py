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
    
    # READ INPUTS
    string1 = 'Usage: python latticeIndexing.py'
    string2 = ' --runNumber <runNumber>'
    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print string1 + string2 
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print string1 + string2 
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

        
    
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
            if (j.selectionFlag == 1 and 
                j.runNumber == runNumber): 
                nProcessedImages = nProcessedImages + 1
                myPath = '%s/latticeDictionary_r%s_image_%s.pkl'%(outputPath, 
                                                                  j.runNumber, 
                                                                  j.imageNumber)
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
        print 'STORING INDEXING LOG FOR RUN %s IN FILE:'%runNumber 
        print '%s/r%s_allLatticesDictionary.txt\n'%(outputPath, 
                                                    runNumber)
        logFile = '%s/r%s_allLatticesDictionary.txt'%(outputPath, 
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
            str1 = '    h    k    observed radius    observed azimuth'
            str2 = '    observed I    radial delta    azimuth delta'
            str3 = '    obs peak n    predicted peak n'
            str4 = '    predicted radius    predicted azimuth\n'
            f.write(str1+str2+str3+str4)
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
        

if __name__ == "__main__":
    mergeFunction(sys.argv[1:])