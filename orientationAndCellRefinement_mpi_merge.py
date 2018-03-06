# -*- coding: utf-8 -*-

import sys
import getopt
import os
import pickle
import joblib

       
    
def merge(myArguments): 
        
    # READ INPUTS
    str1 = 'Usage: python orientationAndCellRefinement.py --runNumber <runNumber>'
    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print str1 
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print str1 
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        
                                
    
    latticeFolder = '/Output_r%s/LatticeIndexing'%runNumber
    refinementFolder = './Output_r%s/OrientationAndCellRefinement'%runNumber       
    
    # LOAD LATTICES
    latticeDictionaryFile = '.%s/r%s_allLatticesDictionary.pkl'%(latticeFolder, 
                                                                 runNumber)      
    if not os.path.exists(latticeDictionaryFile):
        print 'File %s not found.'%latticeDictionaryFile
        return   
        
    fRead = open(latticeDictionaryFile, 'rb')
    myData = pickle.load(fRead)                                                
    fRead.close()       

    # MERGE REFINED LATTICES TO ONE DICTIONARY
    allRefinedLatticesDictionary = {}
    for i, j in myData.items():                                                                  
        myPath = '%s/r%s_Image%s_Lattice%s.jbl'%(refinementFolder, 
                                                 j.runNumber, 
                                                 j.imageNumber, 
                                                 j.latticeNumberInImage)
        if os.path.exists(myPath):                 
            refinedLattice = joblib.load(myPath)
            allRefinedLatticesDictionary['%s_Lattice_%d'
                                          %(refinedLattice.fileName, 
                                            refinedLattice.latticeNumberInImage)] = \
                                            refinedLattice                  
        else:
            print 'File %s not found.'%myPath 
               
    allLatticesDictionaryFile = \
    open('%s/r%s_refinedLatticesDictionary.pkl'%(refinementFolder, 
                                                 runNumber), 
                                                 'wb')
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
    
    # REMOVE INDIVIDUAL LATTICE FILES
    os.system('rm %s/*.jbl'%refinementFolder)
    os.system('rm %s/*.npy'%refinementFolder)

      

if __name__ == "__main__":
    print "\n**** MERGING orientationAndCellRefinement RESULTS ****"
    merge(sys.argv[1:])