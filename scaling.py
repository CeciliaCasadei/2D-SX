# -*- coding: utf-8 -*-
import getopt
import joblib
import random
import time
import numpy
import sys

import scaling_calculateScaleFactor


def scalingFunction(myArguments):
    
    # DEFAULTS
    runNumber = '' 
    nSeeds = 10
    deltaQrodThreshold = 0.003
    n_minThreshold = 8
    nTriangles = 100
    productThreshold = 0.25
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python scaling.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_r%s/transformAndScale'%runNumber
    
    #LOAD LATTICES LIST OF MATRICES: n h k qRod I Icorrected h_transformed k_transformed      
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    nLattices = len(myList)
    
    #SCALE LATTICES WITH RESPECT TO (nSeeds) SEEDS
    LtoSS_vector = []
    n = 0
    for i in range(0, nSeeds):        
        mySeed = random.sample(range(nLattices), 1)
        spotsSeed = myList[mySeed[0]]
        
        if spotsSeed.shape[1] == 8:
            startTime = time.time()
            print 'Good seed choice! (n %d)'%mySeed[0]
            n = n+1
            LtoSi_vector = []
            
            fOpen = open('%s/r%s_scaling_%s.txt'%(outputFolder, runNumber, n), 'w')
            fOpen.write('Total: %s lattices\n'%nLattices)
            fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
            fOpen.write('Min pairs n: %d\n'%n_minThreshold)
            fOpen.write('Product threshold: %f\n'%productThreshold)
            fOpen.write('\nSeed n: %s\n'%mySeed[0])
            
            for firstNeighbor in range(0, nLattices):
                spots1stN = myList[firstNeighbor]  
                if spots1stN.shape[1] == 6:
                    scale = numpy.nan
                    fOpen.write('Lattice %s, Bad 1st N (non oriented)\n'%firstNeighbor)
                    LtoSi_vector.append(scale)
                    continue
                
                print 'New 1st neighbor. Lattice %d'%firstNeighbor                
                nGood = 0
                nBad = 0
                spotsSeed = numpy.asarray(spotsSeed, dtype=numpy.float32)
                spots1stN = numpy.asarray(spots1stN, dtype=numpy.float32)
                n_min, scale_seedTo1stN = scaling_calculateScaleFactor.calculateScaleFactorFunction(spotsSeed, spots1stN, deltaQrodThreshold)
                if n_min < n_minThreshold:
                    scale = numpy.nan
                    fOpen.write('Lattice %s, Bad 1st N\n'%firstNeighbor)
                else:
                    secondShell = random.sample(range(nLattices), nTriangles)
                    for secondNeighbor in secondShell:                    
                        spots2ndN = myList[secondNeighbor]
                        if spots2ndN.shape[1] == 8:
                            spots2ndN = numpy.asarray(spots2ndN, dtype=numpy.float32)
                            n_min, scale_1stNto2ndN = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots1stN, spots2ndN, deltaQrodThreshold)
                            if n_min >= n_minThreshold:
                                
                                
                                n_min, scale_2ndNtoSeed = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots2ndN, spotsSeed, deltaQrodThreshold)
                                if n_min >= n_minThreshold:
                                    product = scale_seedTo1stN*scale_1stNto2ndN*scale_2ndNtoSeed
                                    print product
                                    if abs(product-1) <= productThreshold:
                                        nGood = nGood + 1
                                    else:
                                        nBad = nBad + 1
                                    
                    if nGood+nBad >= 10 and float(nGood)/(nGood+nBad) >= 0.70:
                        scale = float(1)/scale_seedTo1stN
                        fOpen.write('Lattice %s, Scale: %.3f (nGood = %d, nBad = %d)\n'%(firstNeighbor, scale, nGood, nBad))
                    else:
                        scale = numpy.nan
                        fOpen.write('Lattice %s, Scaling: N\A (nGood = %d, nBad = %d)\n'%(firstNeighbor, nGood, nBad)) 
                LtoSi_vector.append(scale) # lattice to seed
                
            LtoSS_vector.append(LtoSi_vector)
            runTime = time.time() - startTime
            fOpen.write('\n It took: %.1f s'%runTime)
            fOpen.close
            
            print len(LtoSS_vector)
            for i in range(0, len(LtoSS_vector)):
                print len(LtoSS_vector[i])
            
    joblib.dump(LtoSS_vector, '%s/r%s-scaling.jbl'%(outputFolder, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING scaling ****"
    scalingFunction(sys.argv[1:])    