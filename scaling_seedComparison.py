# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os

def scaling_seedComparisonFunction(myArguments):
    
    # DEFAULTS:
    outputFolder = ''
    
    # READ INPUTS 
    str_in = '--runNumber <runNumber> --outputFolder <outputFolder>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", 
                                                                 "outputFolder="])
    except getopt.GetoptError:
        print 'Usage: python scaling_seedComparison.py %s'%str_in
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_seedComparison.py %s'%str_in
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--outputFolder":
            outputFolder = value

    if outputFolder == '':
        outputFolder = './Output_r%s/transformAndScale'%runNumber
    
    scalesList = joblib.load('%s/r%s-scaling/r%s-scaling.jbl'
                              %(outputFolder, runNumber, runNumber))  # L to S 
    nSeeds = len(scalesList)
    print 'N seeds: %d'%nSeeds
            
    # FOR EACH SEED, COUNT N OF SCALED LATTICES 
    seedScores = []  
    for i in range(0, nSeeds):
        latticesToSeed = scalesList[i]
        nLattices = len(latticesToSeed)
        nScaled = 0
        for latticeToSeed in latticesToSeed:
            if not numpy.isnan(latticeToSeed):
                nScaled = nScaled + 1
        seedScores.append(nScaled)
    print seedScores
    print 'N lattices: %d'%nLattices   
        
    # ORDER SEEDS FROM BEST TO WORST 
    orderedScalesList = []
    seedScores = numpy.asarray(seedScores)
    print 'Seed:\tN lattices:'
    for i in range(0, nSeeds):
        if len(numpy.argwhere(seedScores > 10)) == 0: # was 20
            break
        bestSeed = numpy.argmax(seedScores)
        print '%d\t%d'%(bestSeed, seedScores[bestSeed])
        seedScores[bestSeed] = 0
        orderedScalesList.append(scalesList[bestSeed])
    nGoodSeeds = len(orderedScalesList)    
    print 'N good seeds: %d'%nGoodSeeds  
    
    # CHECK SEEDS CONSISTENCY (S_i-S_j-L TRIANGLES) AND BUILD S(S_i-S_j) MATRIX  
    S_SeedSeed = numpy.zeros(shape=(nGoodSeeds, nGoodSeeds))    
    for seed1 in range(0, nGoodSeeds):
        for seed2 in range(0, nGoodSeeds):
        #for seed2 in range(seed1+1, nGoodSeeds):
            S1toS2 = []
            for lattice in range(0, nLattices):
                if (not numpy.isnan(orderedScalesList[seed1][lattice]) and 
                    not numpy.isnan(orderedScalesList[seed2][lattice])):
                    latticeToS1 = orderedScalesList[seed1][lattice]
                    latticeToS2 = orderedScalesList[seed2][lattice]
                    S1toS2_item = latticeToS2 / latticeToS1
                    S1toS2.append(S1toS2_item)
            print '\nSeeds: %d %d'%(seed1, seed2)
            N = len(S1toS2)
            print 'Common scaled lattices: %d'%N
            if N < 5:
                S_SeedSeed[seed1, seed2] = numpy.nan
            else:            
                total = 0
                for i in S1toS2:
                    total = total + i
                averageS1toS2 = total / N
                
                sumOfSquares = 0
                for i in S1toS2:
                    sumOfSquares = sumOfSquares + ((i - averageS1toS2)**2)
                standardDeviation = numpy.sqrt(sumOfSquares/N)
                print 'Scale seed %d to seed %d: %.3f +- %.3f'%(seed1, 
                                                                seed2, 
                                                                averageS1toS2, 
                                                                standardDeviation)
                print 'Relative error: %.2f'%(standardDeviation/averageS1toS2)
                S_SeedSeed[seed1, seed2] = averageS1toS2
                    
    numpy.set_printoptions(precision=3)
    print '\nS_SeedSeed:'
    print S_SeedSeed
    print '\nS_SeedSeed * S_SeedSeed.T'  
    print S_SeedSeed*S_SeedSeed.T
                    
              
    
    # CHECK SEEDS CONSISTENCY (Seed_i-Seed_j-Seed_k TRIANGLES)
    products = []
    for i in range(0, nGoodSeeds):
        for j in range(i+1, nGoodSeeds):
            if numpy.isnan(S_SeedSeed[i,j]):
                continue
            for k in range(0, nGoodSeeds):
                if numpy.isnan(S_SeedSeed[j,k]) or numpy.isnan(S_SeedSeed[k,i]):
                    continue
                if k != i:
                    product = S_SeedSeed[i,j]*S_SeedSeed[j,k]*S_SeedSeed[k,i]
                    products.append(product)
    total = 0
    for i in products:
        total = total + i
    if len(products) != 0:
        average = total / len(products)
        sumOfSquares = 0
        for i in products:
            sumOfSquares = sumOfSquares + (i-average)**2
            productsError = numpy.sqrt(sumOfSquares/len(products))
        print '\nSeed - seed - seed triangles: %.3f +- %.3f'%(average, 
                                                              productsError)
    
    
    
    # SCALE ALL LATTICES WITH RESPECT TO REFERENCE SEED [0]                
    latticeScales = []
    for lattice in range(0, nLattices):
        flag = 0
        for seed in range(0, nGoodSeeds):
            if (not numpy.isnan(orderedScalesList[seed][lattice]) and 
                not numpy.isnan(S_SeedSeed[seed, 0])):
                # S(latticeToSeed) * S(seedToReferenceSeed)
                latticeToRefSeed \
                = orderedScalesList[seed][lattice] * S_SeedSeed[seed, 0] 
                latticeScales.append(latticeToRefSeed)
                flag = 1
                break
        if flag == 0:
            latticeScales.append(numpy.nan)
    
    if not os.path.exists('%s/r%s-finalScales'%(outputFolder, runNumber)):
        os.mkdir('%s/r%s-finalScales'%(outputFolder, runNumber))
    joblib.dump(latticeScales, 
                '%s/r%s-finalScales/r%s-finalScales.jbl'%(outputFolder, 
                                                          runNumber, 
                                                          runNumber))
    print 'N SCALE FACTORS BEFORE NORMALIZATION: ', len(latticeScales)


    # SET AVERAGE SCALE TO 1
    nScaled = 0
    scaleSum = 0
    for scale in latticeScales:
        if not numpy.isnan(scale):
            nScaled = nScaled + 1
            scaleSum = scaleSum + scale
    scaleAvg = scaleSum / nScaled
    print 'AVERAGE SCALE: ', scaleAvg
    latticeScales_normalized = []
    for scale in latticeScales:
        if not numpy.isnan(scale):
            scale = scale / scaleAvg
        latticeScales_normalized.append(scale)
    print 'N SCALE FACTORS AFTER NORMALIZATION: ', len(latticeScales_normalized)
    
    # CHECK NORMALIZATION:
    nScaled = 0
    scaleSum = 0
    for scale in latticeScales_normalized:
        if not numpy.isnan(scale):
            nScaled = nScaled + 1
            scaleSum = scaleSum + scale
    scaleAvg = scaleSum / nScaled
    print 'AVERAGE SCALE (AFTER NORMALIZATION): ', scaleAvg
    
    joblib.dump(latticeScales_normalized, 
                '%s/r%s-finalScales/r%s-finalScales-normalized.jbl'%(outputFolder, 
                                                                     runNumber, 
                                                                     runNumber))
    
    
    # LOG FINAL SCALES   
    fOpen = open('%s/r%s-finalScales.txt'%(outputFolder, runNumber), 'w')
    n = 0
    nScaled_final = 0
    for scale in latticeScales_normalized:
        fOpen.write('Lattice: %d - Scale: %s \n'%(n, scale))
        n = n+1
        if not numpy.isnan(scale):
            nScaled_final = nScaled_final + 1
    percentageScaled = float(nScaled_final)/nLattices
    fOpen.write('Fraction of scaled lattices: %d/%d = %.3f'%(nScaled_final,
                                                             nLattices,
                                                             percentageScaled))
    fOpen.close()
    print '\nFraction of scaled lattices: %d/%d = %.3f'%(nScaled_final,
                                                         nLattices,
                                                         percentageScaled)

if __name__ == "__main__":
    print "\n**** CALLING scaling_seedComparison ****"
    scaling_seedComparisonFunction(sys.argv[1:])    