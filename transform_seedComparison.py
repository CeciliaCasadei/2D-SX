# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy

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
    elif index == 3:
        matrix = inversion_permutation
    else:
        matrix = numpy.nan
    return matrix
    
def matrixToIndex(matrix):
    identity = numpy.matrix([[1, 0],[0, 1]])
    inversion = numpy.matrix([[-1, 0],[0, -1]])
    permutation = numpy.matrix([[0, 1],[1, 0]])
    inversion_permutation = numpy.matrix([[0, -1],[-1, 0]])    
    if numpy.array_equal(matrix, identity):
        index = 0
    elif numpy.array_equal(matrix, inversion):
        index = 1
    elif numpy.array_equal(matrix, permutation):
        index = 2
    elif numpy.array_equal(matrix, inversion_permutation):
        index = 3
    else:
        index = -1
    return index
    
def transform_seedComparisonFunction(myArguments):

    runNumber = ''
    
    # READ INPUTS  
    str_in = '--runNumber <runNumber> --nSeeds <nSeeds> --dQrod <dQrod>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", 
                                                                 "nSeeds=", 
                                                                 "dQrod="])
    except getopt.GetoptError:
        print 'Usage: python transform_seedComparison.py %s'%str_in
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_seedComparison.py %s'%str_in
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--nSeeds":
            nSeeds = int(value)
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        

    outputFolder = './Output_r%s/transformAndScale'%runNumber
    orientationsList = joblib.load('%s/r%s-%dseeds-dQrod%.3f-orientations.jbl'
                                    %(outputFolder, 
                                      runNumber, 
                                      nSeeds, 
                                      deltaQrodThreshold))
   
    # FOR EACH SEED, COUNT N OF ORIENTED LATTICES 
    seedScores = []  
    for i in range(0, nSeeds):
        seedToLattices = orientationsList[i]
        nLattices = len(seedToLattices)
        nOriented = 0
        for seedToLattice in seedToLattices:
            if not numpy.isnan(seedToLattice):
                nOriented = nOriented + 1
        seedScores.append(nOriented)
    print 'N oriented lattices:'
    print seedScores   
    
    # ORDER SEEDS FROM BEST TO WORST 
    orderedOrientationsList = []
    seedScores = numpy.asarray(seedScores)
    print 'Seed:\tN lattices:'
    for i in range(0, nSeeds):
        if len(numpy.argwhere(seedScores > 20)) == 0:
            break
        bestSeed = numpy.argmax(seedScores)
        print '%d\t%d'%(bestSeed, seedScores[bestSeed])
        seedScores[bestSeed] = 0
        orderedOrientationsList.append(orientationsList[bestSeed])
          
    nGoodSeeds = len(orderedOrientationsList)    
    print 'N good seeds: %d'%nGoodSeeds
        
    # CHECK SEEDS CONSISTENCY (Seed_i-Seed_j-L TRIANGLES)
    print '\nChecking seeds consistency (I): Seed_i-Seed_j-L TRIANGLES'
    for seed1 in range(0, nGoodSeeds):
        for seed2 in range(seed1+1, nGoodSeeds):
            S1toS2 = []
            for lattice in range(0, nLattices):
                if (not numpy.isnan(orderedOrientationsList[seed1][lattice]) and 
                    not numpy.isnan(orderedOrientationsList[seed2][lattice])):
                    LtoS1_matrix = indexToMatrix(orderedOrientationsList[seed1]
                                                                        [lattice])
                    LtoS2_matrix = indexToMatrix(orderedOrientationsList[seed2]
                                                                        [lattice])
                    S1toS2_matrix = LtoS1_matrix * LtoS2_matrix
                    S1toS2_index = matrixToIndex(S1toS2_matrix)
                    S1toS2.append(S1toS2_index)
            print 'Seeds: %d %d'%(seed1, seed2)
            print S1toS2
            for S1toS2_item in S1toS2:
                if not S1toS2_item == S1toS2[0]:
                    print 'WARNING: CONFLICTING SEEDS'
                
    # BUILD T(Seed_i-Seed_j) MATRIX             
    T_SeedSeed = numpy.zeros(shape=(nGoodSeeds, nGoodSeeds))
    for seed1 in range(0, nGoodSeeds):
        for seed2 in range(0, nGoodSeeds):
            T_SeedSeed[seed1, seed2] = numpy.nan
            for lattice in range(0, nLattices):
                if (not numpy.isnan(orderedOrientationsList[seed1][lattice]) and 
                    not numpy.isnan(orderedOrientationsList[seed2][lattice])):
                    LtoS1_matrix = indexToMatrix(orderedOrientationsList[seed1]
                                                                        [lattice])
                    LtoS2_matrix = indexToMatrix(orderedOrientationsList[seed2]
                                                                        [lattice])
                    S1toS2_matrix = LtoS1_matrix * LtoS2_matrix
                    S1toS2_index = matrixToIndex(S1toS2_matrix)
                    T_SeedSeed[seed1, seed2] = S1toS2_index
                    break
    print '\nT(Seed_i-Seed_j) matrix:\n'
    print T_SeedSeed
    
    # CHECK SEEDS CONSISTENCY (Seed_i-Seed_j-Seed_k TRIANGLES)
    print '\nChecking seeds consistency (I): Seed_i-Seed_j-Seed_k TRIANGLES'
    identity = numpy.matrix([[1, 0],[0, 1]])
    for i in range(0, nGoodSeeds):
        for j in range(i, nGoodSeeds):
            if numpy.isnan(T_SeedSeed[i,j]):
                continue
            T_ij = indexToMatrix(T_SeedSeed[i,j])
            for k in range(0, nGoodSeeds):
                if numpy.isnan(T_SeedSeed[i,k]) or numpy.isnan(T_SeedSeed[j,k]):
                    continue
                T_ik = indexToMatrix(T_SeedSeed[i,k])
                T_jk = indexToMatrix(T_SeedSeed[j,k])
                if not numpy.array_equal(T_ij*T_ik*T_jk, identity):
                    print 'PROBLEM: Seed %d - Seed %d - Seed %d'%(i, j, k)
                
    # ORIENT ALL LATTICES WITH RESPECT TO REFERENCE SEED [0]   
    print '\nOrienting all lattices with respect to the same seed.'
    latticeOrientations = []
    for lattice in range(0, nLattices):
        flag = 0
        for seed in range(0, nGoodSeeds):
            if not numpy.isnan(orderedOrientationsList[seed][lattice]):
                seedToLattice = indexToMatrix(orderedOrientationsList[seed]
                                                                     [lattice])
                seedToRefSeed = indexToMatrix(T_SeedSeed[0, seed])
                latticeToRefSeed = matrixToIndex(seedToLattice*seedToRefSeed)
                latticeOrientations.append(latticeToRefSeed)
                flag = 1
                break
        if flag == 0:
            latticeOrientations.append(numpy.nan)
    joblib.dump(latticeOrientations, '%s/r%s-finalOrientations.jbl'%(outputFolder, 
                                                                     runNumber))
    
    # LOG FINAL ORIENTATIONS        
    fOpen = open('%s/r%s-finalOrientations.txt'%(outputFolder, runNumber), 'w')
    n = 0
    nOriented_final = 0
    for orientation in latticeOrientations:
        fOpen.write('Lattice: %d - Orientation: %s \n'%(n, orientation))
        n = n+1
        if not numpy.isnan(orientation):
            nOriented_final = nOriented_final + 1
    percentageOriented = float(nOriented_final)/nLattices
    fOpen.write('Oriented lattices: %d/%d = %.3f '%(nOriented_final,
                                                    nLattices,
                                                    percentageOriented))
    fOpen.close()
    print '\nFraction of oriented lattices: %d/%d = %.3f'%(nOriented_final,
                                                           nLattices,
                                                           percentageOriented)

if __name__ == "__main__":
    print "\n**** CALLING transform_seedComparison ****"
    transform_seedComparisonFunction(sys.argv[1:])