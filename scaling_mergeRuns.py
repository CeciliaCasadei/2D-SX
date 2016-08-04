# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import getopt
import joblib
import random
import numpy
import sys
import os

import scaling_calculateScaleFactor_runMerging


def scaling_mergeRunsFunction(myArguments):
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']

    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["dQrod=", "nMin=", "nLatticePairs="])
    except getopt.GetoptError:
        print 'Usage: python scaling_mergeRuns.py --dQrod <dQrod> --nMin <nMin> --nLatticePairs <nLatticePairs>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_mergeRuns.py --dQrod <dQrod> --nMin <nMin> --nLatticePairs <nLatticePairs>'
            sys.exit()
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nLatticePairs":
            nLatticePairs = int(value)    
            
    nRuns = len(runNumbers)
    scales_RR = numpy.zeros(shape=(nRuns, nRuns))
    for i in range(0, nRuns):
        run_1 = runNumbers[i]
        folder_1 = './Output_runMerging/spotsMatricesList-Transformed-r%s'%run_1
        myList_1 = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'%(folder_1, run_1)) # n h k qRod Iscaled
        nLattices_1 = len(myList_1)
        for j in range(0, nRuns):
            run_2 = runNumbers[j]
            print '\nRun %s Run %s'%(run_1, run_2)
                
            folder_2 = './Output_runMerging/spotsMatricesList-Transformed-r%s'%run_2        
            myList_2 = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'%(folder_2, run_2)) # n h k qRod Iscaled
            nLattices_2 = len(myList_2)
            print '%d %d'%(nLattices_1, nLattices_2)
            
            # SCALE RUN2 WRT RUN1
            R1toR2_vector = []
            for latticePair in range(0, nLatticePairs):                   
                L1_index = random.sample(range(nLattices_1), 1)
                L2_index = random.sample(range(nLattices_2), 1)
                L1 = myList_1[L1_index[0]]
                L2 = myList_2[L2_index[0]]
                
                L1 = numpy.asarray(L1, dtype=numpy.float32) # n h k qRod Iscaled
                L2 = numpy.asarray(L2, dtype=numpy.float32) # n h k qRod Iscaled
                n_min, scale_L1toL2 = scaling_calculateScaleFactor_runMerging.calculateScaleFactorFunction(L1, L2, deltaQrodThreshold) # I2 = scale * I1
                if n_min >= n_minThreshold:
                    R1toR2_vector.append(scale_L1toL2)
            
            N = len(R1toR2_vector)
            if N == 0:
                scales_RR[i, j] = numpy.nan
            else:            
                total = 0
                for item in R1toR2_vector:
                    total = total + item
                averageR1toR2 = total / N
                
                sumOfSquares = 0
                for item in R1toR2_vector:
                    sumOfSquares = sumOfSquares + ((item - averageR1toR2)**2)
                standardDeviation = numpy.sqrt(sumOfSquares/N)
                print 'Scale run %s to run %s: %.3f +- %.3f'%(run_1, run_2, averageR1toR2, standardDeviation)
                print 'Relative error: %.2f'%(standardDeviation/averageR1toR2)
                scales_RR[i, j] = averageR1toR2
     
    print '\nRun - Run scale factors matrix (scales_RR):'           
    print scales_RR
    print '\nscales_RR * transpose(scales_RR):'  
    print scales_RR * scales_RR.T
    
    # CHECK RUN SCALES CONSISTENCY (Run_i-Run_j-Run_k TRIANGLES)
    products = []
    for i in range(0, nRuns):
        for j in range(i+1, nRuns):
            if numpy.isnan(scales_RR[i,j]):
                continue
            for k in range(0, nRuns):
                if numpy.isnan(scales_RR[j,k]) or numpy.isnan(scales_RR[k,i]):
                    continue
                if k == i or k == j:
                    continue
                product = scales_RR[i,j]*scales_RR[j,k]*scales_RR[k,i]
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
    print '\nRun - run - run triangles: %.3f +- %.3f'%(average, productsError)
    
    
    # APPLY RUN SCALE FACTORS
    runToRunScales = scales_RR[0, :]
    
    ###
    runToRunScales = numpy.asarray(runToRunScales)
    print runToRunScales
    avgScale = numpy.average(runToRunScales)
    print avgScale
    runToRunScales = runToRunScales/avgScale
    print runToRunScales
    ###    
    
    for runIndex in range(0, nRuns):
        runNumber = runNumbers[runIndex]
        runScale = runToRunScales[runIndex] # Scale from run_0 to run_runIndex: L_runIndex = scale * L_0
        if not runScale == 0:
            runScale = float(1)/runScale
            latticesList = joblib.load('./Output_runMerging/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(runNumber, runNumber)) # n h k qRod Iscaled_old
            scaledRun = []
            for lattice in latticesList:
                scaledLattice = []
                for spot in lattice:
                    I_scaled = runScale * spot[4]
                    scaled_spot = [spot[0], spot[1], spot[2], spot[3], spot[4], I_scaled] # n h k qRod Iscaled_old I_scaled
                    scaledLattice.append(scaled_spot)
                scaledLattice = numpy.asarray(scaledLattice)
                scaledRun.append(scaledLattice)
            scaleOutputFolder = './Output_runMerging/spotsMatricesList-Scaled-r%s-AvgTo1'%runNumber
            if not os.path.exists(scaleOutputFolder):
                os.mkdir(scaleOutputFolder)
            joblib.dump(scaledRun, '%s/r%s_scaledSpotsMatricesList.jbl'%(scaleOutputFolder, runNumber))
           
if __name__ == "__main__":
    print "\n**** CALLING scaling_mergeRuns ****"
    scaling_mergeRunsFunction(sys.argv[1:])    