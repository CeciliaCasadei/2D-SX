### INSTRUCTIONS ###
# RUN e.g.:    python transform_CCmethod_main.py -h                                                   (see help)
# RUN e.g.:    python transform_CCmethod_main.py --runNumber 0195 --dQrod 0.005 --dSelf 0 --nMin 6    (set some parameters)
# RUN e.g.:    python transform_CCmethod_main.py                                                      (all default parameters)



# -*- coding: utf-8 -*-
import joblib
import numpy 
import random
import time
import sys
import getopt

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
    
    #DEFAULTS:
    runNumber = '0195'
    deltaQrodThreshold = 0.005
    deltaSelfThreshold = 0
    n_minThreshold = 6
    nSeeds = 1
    nUsedLattices = 2
    nTriangles = 100
    nGoodFraction = 0.7
    
    #READ COMMAND LINE ARGUMENTS
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=", "dQrod=", "dSelf=", "nMin=", "nSeeds=", "nLattices=", "nTriangles=", "nGoodFraction="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_main.py --runNumber <runNumber> --dQrod <dQrod> --dSelf <dSelf> --nMin <nMin> --nSeeds <nSeeds> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_main.py --runNumber <runNumber> --dQrod <dQrod> --dSelf <dSelf> --nMin <nMin> --nSeeds <nSeeds> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--dSelf":
            deltaSelfThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nSeeds":
            nSeeds = int(value)
        elif option == "--nLattices":
            nUsedLattices = value
        elif option == "--nTriangles":
            nTriangles = int(value)
        elif option == "--nGoodFraction":
            nGoodFraction = float(value)
            
    transformationFolder = './Output_r%s/transformAndScale'%runNumber
    
    #LOAD LATTICES LIST OF MATRICES: n h k qRod I Icorrected       
    myList = joblib.load('%s/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'%(transformationFolder, runNumber, runNumber))
    nLattices = len(myList)
    if nUsedLattices == 'all':
        nUsedLattices = int(nLattices)
    else:
        nUsedLattices = int(nUsedLattices)
    print 'n lattices to transform: %d'%nUsedLattices
    
    #ORIENT LATTICES 0 - nUsedLattices-1 WITH RESPECT TO (nSeeds) SEEDS
    LtoSS_vector = []
    n = 0
    for i in range(0, 100):        
        mySeed = random.sample(range(nLattices), 1)
        spotsSeed = myList[mySeed[0]]
        
        n_min, avg_CCs = transform_calculateCCs.determineTransformation(spotsSeed, spotsSeed, deltaQrodThreshold)        
        if (avg_CCs[0] - avg_CCs[1]) > 0.04 and (avg_CCs[0] - avg_CCs[2]) > 0.04 and (avg_CCs[0] - avg_CCs[3]) > 0.04 and n_min >= n_minThreshold:
            startTime = time.time()
            print 'Good seed choice! (n %d)'%mySeed[0]
            n = n+1
            LtoSi_vector = []
            
            fOpen = open('%s/TEST_r%s-lattices0to%d-%dseeds-dSelf%.3f-orientations_%s.txt'%(transformationFolder, runNumber, nUsedLattices-1, nSeeds, deltaSelfThreshold, n), 'w')
            fOpen.write('Total: %s lattices\n'%nLattices)
            fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
            fOpen.write('Delta self: %f\n'%deltaSelfThreshold)
            fOpen.write('Min pairs n: %d\n'%n_minThreshold)
            fOpen.write('\nSeed n: %s\n'%mySeed[0])
            
            for firstNeighbor in range(0, nUsedLattices):
                print 'New 1st neighbor. Lattice %d'%firstNeighbor                
                nGood = 0
                nBad = 0
                
                spots1stN = myList[firstNeighbor]               
                n_min, avg_CCs = transform_calculateCCs.determineTransformation(spotsSeed, spots1stN, deltaQrodThreshold)
                if n_min < n_minThreshold:
                    transformation_SeedTo1stN = numpy.nan
                    fOpen.write('Lattice %s, Bad 1st N\n'%firstNeighbor)
                    print 'Lattice %s, Bad 1st N\n'%firstNeighbor
                else:
                    transformation_SeedTo1stN = avg_CCs.index(max(avg_CCs))
                    transformationMatrix_SeedTo1stN = indexToMatrix(transformation_SeedTo1stN)
                    
                    secondShell = random.sample(range(nLattices), nTriangles)
                    for secondNeighbor in secondShell:
                        
                        spots2ndN = myList[secondNeighbor]
                        n_min, avg_CCs = transform_calculateCCs.determineTransformation(spots2ndN, spots2ndN, deltaQrodThreshold)
                        if (avg_CCs[0] - avg_CCs[1]) > deltaSelfThreshold and (avg_CCs[0] - avg_CCs[2]) > deltaSelfThreshold and (avg_CCs[0] - avg_CCs[3]) > deltaSelfThreshold and n_min >= n_minThreshold:
                            n_min, avg_CCs = transform_calculateCCs.determineTransformation(spots1stN, spots2ndN, deltaQrodThreshold)
                            if n_min >= n_minThreshold:
                                transformation_1stNto2ndN = avg_CCs.index(max(avg_CCs))
                                transformationMatrix_1stNto2ndN = indexToMatrix(transformation_1stNto2ndN)
                                
                                n_min, avg_CCs = transform_calculateCCs.determineTransformation(spotsSeed, spots2ndN, deltaQrodThreshold)
                                if n_min >= n_minThreshold:
                                    transformation_2ndNtoSeed = avg_CCs.index(max(avg_CCs))
                                    transformationMatrix_2ndToSeed = indexToMatrix(transformation_2ndNtoSeed)
                                    
                                    if numpy.array_equal(transformationMatrix_SeedTo1stN*transformationMatrix_1stNto2ndN*transformationMatrix_2ndToSeed, identity):
                                        nGood = nGood + 1
                                    else:
                                        nBad = nBad + 1
                        if nGood+nBad >= 0.4*nTriangles and float(nGood)/(nGood+nBad) >= nGoodFraction:
                            break
                                    
                    if nGood+nBad >= 0.4*nTriangles and float(nGood)/(nGood+nBad) >= nGoodFraction:
                        fOpen.write('Lattice %s, Orientation: %s (good = %d, bad = %d)\n'%(firstNeighbor, transformation_SeedTo1stN, nGood, nBad))
                        print 'Lattice %s, Orientation: %s (good = %d, bad = %d)\n'%(firstNeighbor, transformation_SeedTo1stN, nGood, nBad)
                    else:
                        transformation_SeedTo1stN = numpy.nan
                        fOpen.write('Lattice %s, Orientation: N\A (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)) 
                        print 'Lattice %s, Orientation: N\A (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)
                LtoSi_vector.append(transformation_SeedTo1stN)
                
            LtoSS_vector.append(LtoSi_vector)
            runTime = time.time() - startTime
            fOpen.write('\n It took: %.1f s'%runTime)
            fOpen.close
            
        if len(LtoSS_vector) == nSeeds:
            break
        
    joblib.dump(LtoSS_vector, '%s/TEST_r%s-%dseeds-dQrod%.3f-dSelf%.3f-orientations.jbl'%(transformationFolder, runNumber, nSeeds, deltaQrodThreshold, deltaSelfThreshold))
    
if __name__ == "__main__":
    print "\n**** CALLING transform_CCmethod_main ****"
    main(sys.argv[1:])