# -*- coding: utf-8 -*-
import joblib
import numpy 
import random
import sys
import getopt
from mpi4py import MPI

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
    
    #READ COMMAND LINE ARGUMENTS
    str_1 = '--runNumber <runNumber> --dQrod <dQrod>' 
    str_2 = '--nMin <nMin> --nTriangles <nTriangles>'  
    str_3 = '--nGoodFraction <nGoodFraction> --seed <seed>' 
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=", 
                                                               "dQrod=", 
                                                               "nMin=", 
                                                               "nTriangles=", 
                                                               "nGoodFraction=",
                                                               "seed="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_main.py %s %s %s'%(str_1, 
                                                                   str_2, 
                                                                   str_3)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_main.py %s %s %s'%(str_1, 
                                                                       str_2, 
                                                                       str_3)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nTriangles":
            nTriangles = int(value)
        elif option == "--nGoodFraction":
            nGoodFraction = float(value)
        elif option == "--seed":
            seed = int(value)
            
    transformationFolder = './Output_r%s/transformAndScale'%runNumber
    
    #LOAD LATTICES LIST OF MATRICES: h k qRod Icorrected       
    myList = joblib.load('%s/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'
                          %(transformationFolder, runNumber, runNumber))
    nLattices = len(myList)
    spotsSeed = myList[seed]
    
    # LOOP ON LATTICES, DETERMINE EACH LATTICE TRANSFORMATION WRT SEED
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numberOfCores = comm.Get_size()
    if rank == 0:
        print "I see %d cores are available"%(numberOfCores)
             
    for firstNeighbor in range(0, nLattices):
        
        if (firstNeighbor+1)%numberOfCores != rank:
            continue
                     
        nGood = 0
        nBad  = 0
        flag  = 0
        
        # DETERMINE LATTICE TRANSFORMATION
        spots1stN = myList[firstNeighbor]               
        n_min, \
        avg_CCs = transform_calculateCCs.determineTransformation(spotsSeed, 
                                                                 spots1stN, 
                                                                 deltaQrodThreshold)
        
        # GOOD LATTICE
        if n_min >= n_minThreshold: 
            
            transformation_SeedTo1stN = avg_CCs.index(max(avg_CCs))
            transformationMatrix_SeedTo1stN \
            = indexToMatrix(transformation_SeedTo1stN)
                            
            # VERIFY LATTICE TRANSFORMATION BUILDING 3-LATTICE RELATIONSHIPS
            secondShell = random.sample(range(nLattices), nTriangles)
            for secondNeighbor in secondShell:                   
                spots2ndN = myList[secondNeighbor]                   
                n_min, \
                avg_CCs = \
                transform_calculateCCs.determineTransformation(spots1stN, 
                                                               spots2ndN, 
                                                               deltaQrodThreshold)
                if n_min >= n_minThreshold:
                    transformation_1stNto2ndN = avg_CCs.index(max(avg_CCs))
                    transformationMatrix_1stNto2ndN \
                    = indexToMatrix(transformation_1stNto2ndN)                        
                    n_min, \
                    avg_CCs \
                    = transform_calculateCCs.determineTransformation(spotsSeed, 
                                                                     spots2ndN, 
                                                                     deltaQrodThreshold)
                    if n_min >= n_minThreshold:
                        transformation_2ndNtoSeed = avg_CCs.index(max(avg_CCs))
                        transformationMatrix_2ndToSeed \
                        = indexToMatrix(transformation_2ndNtoSeed)
                        
                        product = (transformationMatrix_SeedTo1stN*
                                   transformationMatrix_1stNto2ndN*
                                   transformationMatrix_2ndToSeed)
                                   
                        if numpy.array_equal(product, identity):
                            nGood = nGood + 1
                        else:
                            nBad = nBad + 1
                            
        nTot = nGood + nBad                    
        if (float(nTot)/nTriangles > 0.6 and 
            float(nGood)/nTot > nGoodFraction):
                
            flag = 1

        if flag == 1:
            result = [firstNeighbor, transformation_SeedTo1stN]
        else:
            result = [firstNeighbor, numpy.nan]
          
        joblib.dump(result, 
                    '%s/r%s_orientations_L_%s.jbl'
                     %(transformationFolder, 
                       runNumber,
                       str(firstNeighbor).zfill(5)))
    
if __name__ == "__main__":
    #print "\n**** CALLING transform_CCmethod_main ****"
    main(sys.argv[1:])