# -*- coding: utf-8 -*-
import getopt
import joblib
import random
import time
import numpy
import sys
import os

import scaling_calculateScaleFactor


def scalingFunction(myArguments):
    
    # FIXED
    nSeeds = 6
    nTriangles = 100
    
    # DEFAULTS    
    deltaQrodThreshold = 0.003
    productThreshold = 0.25
    resolution_3D = 6.5
    n_minThreshold = 8
        
    input_str_1 = '--runNumber <runNumber>  --dQrod <dQrod>'
    input_str_2 = '--outputFolder <outputFolder> --productThreshold <productThreshold>'
    input_str_3 = '--resolution_3D <resolution_3D> --n_minThreshold <n_minThreshold>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", 
                                             ["runNumber=", 
                                              "dQrod=", 
                                              "outputFolder=", 
                                              "productThreshold=", 
                                              "resolution_3D=",
                                              "n_minThreshold="])
    except getopt.GetoptError:
        print 'Usage: python scaling.py %s %s %s'%(input_str_1, input_str_2, input_str_3)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling.py %s %s %s'%(input_str_1, input_str_2, input_str_3)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--outputFolder":
            outputFolder = value
        elif option == "--productThreshold":
            productThreshold = float(value)
        elif option == "--resolution_3D":
            resolution_3D = float(value)
        elif option == "--n_minThreshold":
            n_minThreshold = int(value)
    
    if not 'outputFolder' in locals():                    
        outputFolder = './Output_r%s/transformAndScale'%runNumber
        
    print 'OUTPUT FOLDER: ', outputFolder
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod I flag    
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'
                         %(outputFolder, runNumber, runNumber))
    nLattices = len(myList)
    
    # SCALE LATTICES WITH RESPECT TO (nSeeds) SEEDS
    LtoSS_vector = []
    n = 0
    for i in range(0, 100):        
        mySeed = random.sample(range(nLattices), 1)
        spotsSeed = myList[mySeed[0]]
        if spotsSeed[0, 4] == 1:  # check flag value
            n = n + 1
            startTime = time.time()
            
            print 'Good seed choice! (n %d)'%mySeed[0]          
            fOpen = open('%s/r%s_scaling_%s.txt'%(outputFolder, runNumber, n), 'w')
            fOpen.write('Total: %s lattices\n'%nLattices)
            fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
            fOpen.write('Min pairs n: %d\n'%n_minThreshold)
            fOpen.write('Product threshold: %f\n'%productThreshold)
            fOpen.write('\nSeed n: %s\n'%mySeed[0])
            
            # LOOP ON LATTICES, FOR EACH LATTICE DETERMINE SCALE WRT SEED
            LtoSi_vector = []
            for firstNeighbor in range(0, nLattices):
                spots1stN = myList[firstNeighbor]  
                
                # BAD LATTICE
                if spots1stN[0, 4] == 0:
                    scale = numpy.nan
                    fOpen.write('Lattice %s: non oriented.\n'%firstNeighbor)
                    print 'Lattice %s: non oriented.\n'%firstNeighbor
                    
                # GOOD LATTICE
                else:
                    nGood = 0
                    nBad = 0
                    n_min, \
                    scale_seedTo1stN, \
                    I1, \
                    I2 = scaling_calculateScaleFactor.calculateScaleFactorFunction(spotsSeed, 
                                                                                   spots1stN, 
                                                                                   deltaQrodThreshold, 
                                                                                   resolution_3D)
                    # BAD LATTICE
                    if n_min < n_minThreshold:
                        scale = numpy.nan
                        fOpen.write('Lattice %s, n_min below threshold.\n'%firstNeighbor)
                        print 'Lattice %s, n_min below threshold.\n'%firstNeighbor
                    # VERIFY SCALE
                    else:
                        secondShell = random.sample(range(nLattices), nTriangles)
                        for secondNeighbor in secondShell:                    
                            spots2ndN = myList[secondNeighbor]
                            if spots2ndN[0, 4] == 1:                         
                                n_min, \
                                scale_1stNto2ndN, \
                                I1, \
                                I2 = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots1stN, 
                                                                                               spots2ndN, 
                                                                                               deltaQrodThreshold, 
                                                                                               resolution_3D)
                                if n_min >= n_minThreshold:
                                    n_min, \
                                    scale_2ndNtoSeed, \
                                    I1, \
                                    I2 = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots2ndN, 
                                                                                                   spotsSeed, 
                                                                                                   deltaQrodThreshold,
                                                                                                   resolution_3D)
                                    if n_min >= n_minThreshold:
                                        product = scale_seedTo1stN*scale_1stNto2ndN*scale_2ndNtoSeed
                                        if abs(product-1) <= productThreshold:
                                            nGood = nGood + 1
                                        else:
                                            nBad = nBad + 1
                                        
                        if nGood+nBad >= 10 and float(nGood)/(nGood+nBad) >= 0.70:
                            scale = float(1)/scale_seedTo1stN
                            fOpen.write('Lattice %s, Scale: %.3f (good = %d, bad = %d)\n'%(firstNeighbor, scale, nGood, nBad))
                            print 'Lattice %s, Scale: %.3f (good = %d, bad = %d)\n'%(firstNeighbor, scale, nGood, nBad)
                        else:
                            scale = numpy.nan
                            fOpen.write('Lattice %s, Scale: n/a (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)) 
                            print 'Lattice %s, Scale: n/a (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)
                LtoSi_vector.append(scale) # lattice to seed
                
            LtoSS_vector.append(LtoSi_vector)
            runTime = time.time() - startTime
            fOpen.write('\nIt took: %.1f s'%runTime)
            fOpen.close
        
        if n == nSeeds:
            break
                
    if not os.path.exists('%s/r%s-scaling'%(outputFolder, runNumber)):
        os.mkdir('%s/r%s-scaling'%(outputFolder, runNumber))
        
    joblib.dump(LtoSS_vector, '%s/r%s-scaling/r%s-scaling.jbl'%(outputFolder, runNumber, runNumber))

if __name__ == "__main__":
    print "\n**** CALLING scaling ****"
    scalingFunction(sys.argv[1:])    