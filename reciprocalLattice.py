# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:39:04 2016
@author: casadei_c
CALCULATE RECIPROCAL LATTICE IN REFERENCE ORIENTATION.
"""
import sys
import getopt
import os
import pickle
import buildReciprocalLattice
import plotReciprocalLattice

def reciprocalLatticeFunction(myArguments):
    
    # DEFAULTS
    referenceCellSize = 62.45
    hmax = 100
    kmax = 100
    resolutionLimit = 5.0
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["referenceCellSize=", "hmax=", "kmax=", "resolutionLimit=", "outFolder="])
    except getopt.GetoptError:
        print 'Usage: python reciprocalLattice.py --referenceCellSize <referenceCellSize> --hmax <hmax> --kmax <kmax> --resolutionLimit <resolutionLimit> --outFolder <outFolder>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python reciprocalLattice.py --referenceCellSize <referenceCellSize> --hmax <hmax> --kmax <kmax> --resolutionLimit <resolutionLimit> --outFolder <outFolder>'
            sys.exit()
        elif option == "--referenceCellSize":
            referenceCellSize = float(value)
        elif option == "--hmax":
            hmax = int(value)
        elif option == "--kmax":
            kmax = int(value)
        elif option == "--resolutionLimit":
            resolutionLimit = float(value) 
        elif option == "--outFolder":
            outFolder = value
        
    if not os.path.exists('%s/ReferenceReciprocalLattice'%outFolder):   
        os.mkdir('%s/ReferenceReciprocalLattice'%outFolder)  
    
    # CALCULATE RECIPROCAL LATTICE, SAVE BINARY AND FIGURE
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(referenceCellSize, hmax, kmax, resolutionLimit)
    
    f = open('%s/ReferenceReciprocalLattice/reciprocalLattice_cellSize_%.3f.pkl'%(outFolder, referenceCellSize), 'wb') 
    pickle.dump(reciprocalLattice, f)
    f.close()
    
    plotReciprocalLattice.plotReciprocalLatticeFunction(referenceCellSize, outFolder)
    
    
    
if __name__ == "__main__":
    print "\n**** CALLING reciprocalLattice ****"
    reciprocalLatticeFunction(sys.argv[1:])