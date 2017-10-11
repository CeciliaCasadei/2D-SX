# -*- coding: utf-8 -*-
import getopt
import joblib
import random
import numpy
import sys
import os

import scaling_calculateScaleFactor

def scaling_mergeRunsFunction(myArguments):
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    resolution_3D = 6.5 

    str_input = '--dQrod <dQrod> --nMin <nMin> --nLatticePairs <nLatticePairs> --resolution_3D <resolution_3D>'
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["dQrod=", "nMin=", "nLatticePairs=", "resolution_3D="])
    except getopt.GetoptError:
        print 'Usage: python scaling_mergeRuns.py %s'%str_input
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_mergeRuns.py %s'%str_input
            sys.exit()
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nLatticePairs":
            nLatticePairs = int(value)  
        elif option == "--resolution_3D":
            resolution_3D = float(value)
            
    nRuns = len(runNumbers)
    scales_RR = numpy.zeros(shape=(nRuns, nRuns))
    for i in range(0, nRuns):
        run_1 = runNumbers[i]
        folder_1 = './Output_runMerging/spotsMatricesList-Transformed-r%s'%run_1
        myList_1 = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'%(folder_1, run_1))    # h k qRod I flag i_unassembled j_unassembled
        nLattices_1 = len(myList_1)
        for j in range(0, nRuns):
            run_2 = runNumbers[j]
            print '\nRun %s Run %s'%(run_1, run_2)
                
            folder_2 = './Output_runMerging/spotsMatricesList-Transformed-r%s'%run_2        
            myList_2 = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'%(folder_2, run_2)) # h k qRod I flag i_unassembled j_unassembled
            nLattices_2 = len(myList_2)
            print '%d %d'%(nLattices_1, nLattices_2)
            
            # SCALE RUN2 WRT RUN1
            R1toR2_vector = []
            for latticePair in range(0, nLatticePairs):                   
                L1_index = random.sample(range(nLattices_1), 1)
                L2_index = random.sample(range(nLattices_2), 1)
                L1 = myList_1[L1_index[0]]
                L2 = myList_2[L2_index[0]]
                L1 = numpy.asarray(L1, dtype=numpy.float32) # h k qRod I flag i_unassembled j_unassembled
                L2 = numpy.asarray(L2, dtype=numpy.float32) # h k qRod I flag i_unassembled j_unassembled
                if L1[0, 4] == 1 and L2[0, 4] == 1:
                    n_min, scale_L1toL2, I1, I2 = scaling_calculateScaleFactor.calculateScaleFactorFunction(L1, 
                                                                                                            L2, 
                                                                                                            deltaQrodThreshold, 
                                                                                                            resolution_3D) # I2 = scale * I1
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
    
    
    # EXTRACT RUN SCALE FACTORS
    runToRunScales = scales_RR[:, 1] # Scale run_i to run_1, was run_0 !!!
    
    # NORMALIZE SCALES
    runToRunScales = numpy.asarray(runToRunScales)
    print runToRunScales
    clean_runToRunScales = [runToRunScales[i] for i in range(0, len(runToRunScales)) if not numpy.isnan(runToRunScales[i])]
    avgScale = numpy.average(clean_runToRunScales)
    print avgScale
    runToRunScales = runToRunScales/avgScale
    print runToRunScales
      
    # APPLY RUN SCALE FACTORS
    for runIndex in range(0, nRuns):
        runNumber = runNumbers[runIndex]
        runScale = runToRunScales[runIndex] # Scale from run_runIndex to run_0: L_0 = scale * L_runIndex
        latticesList = joblib.load('./Output_runMerging/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(runNumber, runNumber)) # h k qRod I flag
        scaledRun = []
        for lattice in latticesList:
            lattice = numpy.asarray(lattice)
            if lattice[0, 4] == 0 or numpy.isnan(runScale):
                flag = 0
                scaledLattice = []
                for spot in lattice:
                    scaled_spot = [spot[0], spot[1], spot[2], spot[3], flag, spot[5], spot[6]]    # h k qRod I flag i_unassembled j_unassembled
                    scaledLattice.append(scaled_spot)
            else: 
                flag = 1
                scaledLattice = []
                for spot in lattice:
                    I_scaled = runScale * spot[3]
                    scaled_spot = [spot[0], spot[1], spot[2], I_scaled, flag, spot[5], spot[6]]   # h k qRod I flag i_unassembled j_unassembled
                    scaledLattice.append(scaled_spot)
            scaledLattice = numpy.asarray(scaledLattice)
            scaledRun.append(scaledLattice)
        scaleOutputFolder = './Output_runMerging/spotsMatricesList-Scaled-r%s'%runNumber
        if not os.path.exists(scaleOutputFolder):
            os.mkdir(scaleOutputFolder)
        joblib.dump(scaledRun, '%s/r%s_scaledSpotsMatricesList.jbl'%(scaleOutputFolder, runNumber))
           
if __name__ == "__main__":
    print "\n**** CALLING scaling_mergeRuns ****"
    scaling_mergeRunsFunction(sys.argv[1:])    