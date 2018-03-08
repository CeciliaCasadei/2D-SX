# -*- coding: utf-8 -*-
import sys
import getopt
import numpy
import joblib
import os

def transform_makeSpotsMatrixFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python transform_makeSpotsMatrix.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_makeSpotsMatrix.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_r%s/transformAndScale'%runNumber
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    spotMatricesFolder = '%s/spotsMatricesList-r%s'%(outputFolder, runNumber)
    if not os.path.exists(spotMatricesFolder):
        os.mkdir(spotMatricesFolder)
        
    os.system('ls ./Output_r%s/UnassembledImageProcessing/Ordered*.jbl > %s/list-r%s.txt'%(runNumber, 
                                                                                           spotMatricesFolder, 
                                                                                           runNumber))
    f = open('%s/list-r%s.txt'%(spotMatricesFolder, runNumber), 'r')
    myList = list(f)
    f.close()
    nLattices = len(myList)
    print 'N processed lattices: %d'%nLattices
    
    spotsMatricesList = []
    n = 0
    for myLatticeFile in myList:
        n = n+1
        print n
        myLatticeFile = myLatticeFile.strip('\n')   ### e.g. OrderedIntegratedIntensities_r0195_Img0127_Lattice1.jbl
        spotsMatrix2D = joblib.load(myLatticeFile)  ### h k qRod Icorrected flag=1 i_unassembled j_unassembled
        spotsMatrix2D = numpy.asarray(spotsMatrix2D, dtype=numpy.float32)
        print spotsMatrix2D.shape
        spotsMatricesList.append(spotsMatrix2D)     ### spotsMatricesList[lattice n in list.txt][n spot]
               
    nLattices = len(spotsMatricesList)
    print 'N processed lattices: %d'%nLattices
    
    joblib.dump(spotsMatricesList, '%s/r%s_spotsMatricesList.jbl'%(spotMatricesFolder, 
                                                                   runNumber))

if __name__ == "__main__":
    print "\n**** CALLING transform_makeSpotsMatrix ****"
    transform_makeSpotsMatrixFunction(sys.argv[1:])    