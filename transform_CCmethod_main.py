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
    runNumber = ''
    deltaQrodThreshold = 0.005
    n_minThreshold = 6
    nSeeds = 6
    nUsedLattices = 'all'
    nTriangles = 100
    nGoodFraction = 0.7
    
    #READ COMMAND LINE ARGUMENTS
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=", "dQrod=", "nMin=", "nSeeds=", "nLattices=", "nTriangles=", "nGoodFraction="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_main.py --runNumber <runNumber> --dQrod <dQrod> --nMin <nMin> --nSeeds <nSeeds> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_main.py --runNumber <runNumber> --dQrod <dQrod> --nMin <nMin> --nSeeds <nSeeds> --nLattices <nLattices> --nTriangles <nTriangles> --nGoodFraction <nGoodFraction>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
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
    
    #LOAD LATTICES LIST OF MATRICES: h k qRod Icorrected       
    myList = joblib.load('%s/spotsMatricesList-r%s/r%s_spotsMatricesList.jbl'%(transformationFolder, runNumber, runNumber))
    nLattices = len(myList)
    if nUsedLattices == 'all':
        nUsedLattices = int(nLattices)
    else:
        nUsedLattices = int(nUsedLattices)
    print 'n lattices to transform: %d'%nUsedLattices
    
    #ORIENT LATTICES 0 - nUsedLattices-1 WITH RESPECT TO (nSeeds) SEEDS
    LtoSS_vector = []    
    for i in range(0, nSeeds):
        LtoSi_vector = []
        mySeed = random.sample(range(nLattices), 1)
        spotsSeed = myList[mySeed[0]]
        
        startTime = time.time()
        print 'New seed: %d'%mySeed[0]
        fOpen = open('%s/r%s-%dseeds-orientations_%s.txt'%(transformationFolder, runNumber, nSeeds, i), 'w')
        fOpen.write('Total: %s lattices\n'%nLattices)
        fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)        
        fOpen.write('Min pairs n: %d\n'%n_minThreshold)
        fOpen.write('\nSeed n: %s\n'%mySeed[0])
        
        # LOOP ON LATTICES, DETERMINE EACH LATTICE TRANSFORMATION WRT SEED
        for firstNeighbor in range(0, nUsedLattices):
            print 'Lattice %d'%firstNeighbor                
            nGood = 0
            nBad = 0
            
            # DETERMINE LATTICE TRANSFORMATION
            spots1stN = myList[firstNeighbor]               
            n_min, avg_CCs = transform_calculateCCs.determineTransformation(spotsSeed, spots1stN, deltaQrodThreshold)
            
            # BAD LATTICE
            if n_min < n_minThreshold:
                transformation_SeedTo1stN = numpy.nan
                fOpen.write('Lattice %s: n_min lattice-seed below threshold.\n'%firstNeighbor)
                print 'Lattice %s: n_min lattice-seed below threshold.\n'%firstNeighbor
                
            # GOOD LATTICE
            else:
                transformation_SeedTo1stN = avg_CCs.index(max(avg_CCs))
                transformationMatrix_SeedTo1stN = indexToMatrix(transformation_SeedTo1stN)
                                
                # VERIFY LATTICE TRANSFORMATION BUILDING TRIANGULAR RELATIONSHIPS
                secondShell = random.sample(range(nLattices), nTriangles)
                for secondNeighbor in secondShell:                   
                    spots2ndN = myList[secondNeighbor]                   
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
                    fOpen.write('Lattice %s, Orientation: n/a (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)) 
                    print 'Lattice %s, Orientation: n/a (good = %d, bad = %d)\n'%(firstNeighbor, nGood, nBad)
            LtoSi_vector.append(transformation_SeedTo1stN)
            
        LtoSS_vector.append(LtoSi_vector)
        runTime = time.time() - startTime
        fOpen.write('\nIt took: %.1f s'%runTime)
        fOpen.close        
    joblib.dump(LtoSS_vector, '%s/r%s-%dseeds-dQrod%.3f-orientations.jbl'%(transformationFolder, runNumber, nSeeds, deltaQrodThreshold))
    
if __name__ == "__main__":
    print "\n**** CALLING transform_CCmethod_main ****"
    main(sys.argv[1:])