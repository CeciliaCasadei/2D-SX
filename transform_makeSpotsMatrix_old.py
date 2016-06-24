# -*- coding: utf-8 -*-

import numpy
import joblib
import os

runNumber = '0195'
outputFolder = './Output/transformAndScale'
if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)
spotMatricesFolder = './Output/transformAndScale/spotsMatricesList-r%s'%runNumber
if not os.path.exists(spotMatricesFolder):
    os.mkdir(spotMatricesFolder)
    
os.system('ls ./Output/UnassembledImageProcessing/*.jbl > %s/list-r%s.txt'%(spotMatricesFolder, runNumber))
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
    myLatticeFile = myLatticeFile.strip('\n')  ### e.g. SingleLattice_r0195_Img127_Lattice1.jbl
    myLattice = joblib.load(myLatticeFile)
    spotsMatrix2D = myLattice.orderedIntegratedIntensities
    #spotsMatrix2D = joblib.load(myLatticeFile)  ## e.g. OrderedIntegratedIntensities_r0195_Img127_Lattice1.jbl
    spotsMatrix2D = numpy.asarray(spotsMatrix2D, dtype=numpy.float32)
    print spotsMatrix2D.shape
    spotsMatricesList.append(spotsMatrix2D)
    
    
nLattices = len(spotsMatricesList)
print nLattices
print spotsMatricesList[5][17] #spotsMatrix3D[lattice n in list.txt][n spot]

joblib.dump(spotsMatricesList, '%s/r%s_spotsMatricesList.jbl'%(spotMatricesFolder, runNumber))