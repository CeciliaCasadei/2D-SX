# -*- coding: utf-8 -*-
import joblib
import numpy 
import random
import time
import sys
import getopt
import os

### CYTHON MODULES ###
import transform_calculateCCs

def indexToMatrix(index):
    identity = numpy.matrix([[1, 0],[0, 1]])
    inversion = numpy.matrix([[-1, 0],[0, -1]])
    permutation = numpy.matrix([[0, 1],[1, 0]])
    inversion_permutation = numpy.matrix([[0, -1],[-1, 0]])    
    if index == 0:
        matrix = identity
    elif index == 1:
        matrix = inversion
    elif index == 2:
        matrix = permutation
    else:
        matrix = inversion_permutation
    return matrix
    
def main(myArguments): 
    identity = numpy.matrix([[1, 0],[0, 1]])
    
    # DEFAULTS:
    runNumber = ''
    deltaQrodThreshold = 0.001
    n_minThreshold = 6
    nUsedLattices = 30
    nTriangles = 100
    nGoodFraction = 0.7
    
    # READ COMMAND LINE ARGUMENTS
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=", "dQrod=", "nMin=", "nLattices=", "nTriangles=", "nGoodFraction="])
    except getopt.GetoptError:
        print 'Usage: python model_transformVsModel.py --runNumber <runNumber> --dQrod <dQrod> --nMin <nMin> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_transformVsModel.py --runNumber <runNumber> --dQrod <dQrod> --nMin <nMin> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nLattices":
            nUsedLattices = value
        elif option == "--nTriangles":
            nTriangles = int(value)  
        elif option == "--nGoodFraction":
            nGoodFraction = float(value)
    
    newFolder = './Output_runMergingVsModel'
    #if not os.path.exists(newFolder):
        #os.mkdir(newFolder)         
    transformationFolder = '%s/transformAndScaleToModel_r%s'%(newFolder, runNumber)
    if not os.path.exists(transformationFolder):
        os.mkdir(transformationFolder)

    # LOAD MODEL: h k qRod I
    lattice_model = joblib.load('./Output_runMerging/model/lattice_model.jbl')  
            
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag=1 i_unassembled j_unassembled   
    myList = joblib.load('Output_r%s/transformAndScale/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'%(runNumber, runNumber, runNumber))
    nLattices = len(myList)
    if nUsedLattices == 'all':
        nUsedLattices = int(nLattices)
    else:
        nUsedLattices = int(nUsedLattices)
    print 'Number of lattices to transform: %d'%nUsedLattices
    startTime = time.time()
    
    # ORIENT LATTICES 0 - nUsedLattices-1 WITH RESPECT TO MODEL
    nOriented = 0
    LtoModel_vector = []
    
    fOpen = open('%s/r%s-orientationsVsModel.txt'%(transformationFolder, runNumber), 'w')
    fOpen.write('Total: %s lattices\n'%nLattices)
    fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
    fOpen.write('Min pairs n: %d\n'%n_minThreshold)
   
    for firstNeighbor in range(0, nUsedLattices):
        print '\nLattice %d'%firstNeighbor                
        nGood = 0
        nBad = 0
        
        spots1stN = myList[firstNeighbor]   
        print spots1stN.shape            
        n_min, avg_CCs = transform_calculateCCs.determineTransformation(lattice_model, spots1stN, deltaQrodThreshold)
        if n_min < n_minThreshold:
            transformation_ModelTo1stN = numpy.nan
            fOpen.write('Lattice %s, Bad 1st N\n'%firstNeighbor)
            print 'n_min below threshold (model to 1st neighbour n %s)'%firstNeighbor
        else:
            print 'n_min above threshold (model to 1st neighbour): %d'%n_min
            transformation_ModelTo1stN = avg_CCs.index(max(avg_CCs))
            transformationMatrix_ModelTo1stN = indexToMatrix(transformation_ModelTo1stN)
            
            secondShell = random.sample(range(nLattices), nTriangles)
            for secondNeighbor in secondShell:
                
                spots2ndN = myList[secondNeighbor]
                
                n_min, avg_CCs = transform_calculateCCs.determineTransformation(spots1stN, spots2ndN, 5*deltaQrodThreshold)
                if n_min < n_minThreshold:
                    print 'n_min below threshold (1st to 2nd neighbour)'
                else:
                    transformation_1stNto2ndN = avg_CCs.index(max(avg_CCs))
                    transformationMatrix_1stNto2ndN = indexToMatrix(transformation_1stNto2ndN)
                    
                    n_min, avg_CCs = transform_calculateCCs.determineTransformation(lattice_model, spots2ndN, deltaQrodThreshold)
                    if n_min < n_minThreshold:
                        print 'n_min below threshold (model to 2nd neighbour)'
                    else:
                        transformation_2ndNtoModel = avg_CCs.index(max(avg_CCs))
                        transformationMatrix_2ndToModel = indexToMatrix(transformation_2ndNtoModel)
                        
                        if numpy.array_equal(transformationMatrix_ModelTo1stN*transformationMatrix_1stNto2ndN*transformationMatrix_2ndToModel, identity):
                            nGood = nGood + 1
                        else:
                            nBad = nBad + 1
                            
                if nGood+nBad >= 0.4*nTriangles and float(nGood)/(nGood+nBad) >= nGoodFraction:
                    break
                            
            if nGood+nBad >= 0.4*nTriangles and float(nGood)/(nGood+nBad) >= nGoodFraction:
                fOpen.write('Lattice %s, Orientation: %s (good = %d, bad = %d)\n'%(firstNeighbor, transformation_ModelTo1stN, nGood, nBad))
                print 'Lattice %s, Orientation: %s (good = %d, bad = %d)\n'%(firstNeighbor, transformation_ModelTo1stN, nGood, nBad)
                nOriented = nOriented + 1
            else:
                transformation_ModelTo1stN = numpy.nan
                fOpen.write('Lattice %s, Orientation: N\A (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)) 
                print 'Lattice %s, Orientation: N\A (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)
        LtoModel_vector.append(transformation_ModelTo1stN)
        
    fractionOriented = float(nOriented)/nLattices
    fOpen.write('\nFraction of transformed lattices: %.2f.'%fractionOriented)
    runTime = time.time() - startTime
    fOpen.write('\nIt took: %.1f s'%runTime)
    fOpen.close
            
    joblib.dump(LtoModel_vector, '%s/r%s_orientationsVsModel.jbl'%(transformationFolder, runNumber))
    
if __name__ == "__main__":
    print "\n**** CALLING model_transformVsModel ****"
    main(sys.argv[1:])