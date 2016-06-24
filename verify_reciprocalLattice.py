# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:41:39 2016
@author: casadei_c
PRINT ON TERMINAL RECIPROCAL LATTICE CONTENTS.
"""

import sys
import getopt
import os
import pickle
    
def verifyReciprocalLatticeFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["reciprocalLatticeFile="])
    except getopt.GetoptError:
        print 'Usage: python verify_reciprocalLattice.py --reciprocalLatticeFile <reciprocalLatticeFile>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python verify_reciprocalLattice.py --reciprocalLatticeFile <reciprocalLatticeFile>'
            sys.exit()
        elif option == "--reciprocalLatticeFile":
            reciprocalLatticeFile = value
        
    if not os.path.exists(reciprocalLatticeFile):
        print 'File %s not found.'%reciprocalLatticeFile       
    else:    
        fRead = open(reciprocalLatticeFile, 'rb')
        myData = pickle.load(fRead)
        fRead.close()
        print '\nEXTRACTING DATA FROM %s\n'%reciprocalLatticeFile
        print '    h     k         qx         qy   resolution        q'
        for i in myData:
            print '%5d %5d %10.4f %10.4f %12.2f %8.3f'%(i[0], i[1], i[2], i[3], i[4], i[5])
            
            
        
if __name__ == "__main__":
    print "\n**** VERIFY RECIPROCAL LATTICE ****"
    verifyReciprocalLatticeFunction(sys.argv[1:])