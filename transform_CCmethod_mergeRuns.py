# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
import joblib
import numpy 
import random
import sys
import getopt
import os

### CYTHON MODULES ###
import transform_calculateCCs_mergeRuns



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
    
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    #READ COMMAND LINE ARGUMENTS
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["dQrod=", "nMin=", "nLatticePairs="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_mergeRuns.py --dQrod <dQrod> --nMin <nMin> --nLatticePairs <nLatticePairs>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_mergeRuns.py --dQrod <dQrod> --nMin <nMin> --nLatticePairs <nLatticePairs>'
            sys.exit()
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nLatticePairs":
            nLatticePairs = int(value)
            
    nRuns = len(runNumbers)
    transformations_RR = numpy.zeros(shape=(nRuns, nRuns))
    for i in range(0, nRuns):
        run_1 = runNumbers[i]
        transformationFolder_1 = './Output_r%s/transformAndScale'%run_1
        myList_1 = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(transformationFolder_1, run_1, run_1)) # n h k qRod I Icorrected h_transformed k_transformed Iscaled
        nLattices_1 = len(myList_1)
        for j in range(0, nRuns):
            run_2 = runNumbers[j]
            print '\nRun %s Run %s'%(run_1, run_2)
                
            transformationFolder_2 = './Output_r%s/transformAndScale'%run_2            
            myList_2 = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(transformationFolder_2, run_2, run_2)) # n h k qRod I Icorrected h_transformed k_transformed Iscaled
            nLattices_2 = len(myList_2)
            print '%d %d'%(nLattices_1, nLattices_2)
    
            # ORIENT RUN2 WRT RUN1
            R2toR1_vector = []
            for latticePair in range(0, nLatticePairs):                   
                L1_index = random.sample(range(nLattices_1), 1)
                L2_index = random.sample(range(nLattices_2), 1)
                L1 = myList_1[L1_index[0]]
                L2 = myList_2[L2_index[0]]
                if L1.shape[1] == 9 and L2.shape[1] == 9:
                    L1 = numpy.asarray(L1, dtype=numpy.float32)
                    L2 = numpy.asarray(L2, dtype=numpy.float32)
                    n_min, avg_CCs = transform_calculateCCs_mergeRuns.determineTransformation(L1, L2, deltaQrodThreshold)        
                    if n_min >= n_minThreshold:
                        transformation = avg_CCs.index(max(avg_CCs))
                        R2toR1_vector.append(transformation)
                
            print R2toR1_vector
            n_I = 0
            n_i = 0
            n_p = 0
            n_ip = 0
            for item in R2toR1_vector:
                if item == 0:
                    n_I = n_I + 1
                elif item == 1:
                    n_i = n_i + 1
                elif item == 2:
                    n_p = n_p + 1
                elif item == 3:
                    n_ip = n_ip + 1
                else:
                    print 'ERROR'
            LtoL_transformations = [n_I, n_i, n_p, n_ip]
            transformation_R1_R2 = LtoL_transformations.index(max(LtoL_transformations))
            print transformation_R1_R2
            transformations_RR[i, j] = transformation_R1_R2
    
    print 'Run - run transformations'        
    print transformations_RR
    
    if numpy.array_equal(transformations_RR, transformations_RR.T):
        print 'SYMMETRIC'
        print 'NON SYMMETRIC'
        
    # CHECK RUN TRANSFORMATIONS CONSISTENCY (Run_i-Run_j-Run_k TRIANGLES)
    print '\nChecking run transformations consistency: Run_i-Run_j-Run_k TRIANGLES'
    identity = numpy.matrix([[1, 0],[0, 1]])
    for i in range(0, nRuns):
        for j in range(i, nRuns):
            if numpy.isnan(transformations_RR[i,j]):
                continue
            T_ij = indexToMatrix(transformations_RR[i,j])
            for k in range(0, nRuns):
                if numpy.isnan(transformations_RR[i,k]) or numpy.isnan(transformations_RR[j,k]):
                    continue
                T_ik = indexToMatrix(transformations_RR[i,k])
                T_jk = indexToMatrix(transformations_RR[j,k])
                if not numpy.array_equal(T_ij*T_ik*T_jk, identity):
                    print 'PROBLEM: Run %d - Run %d - Run %d'%(i, j, k)
                else:
                    print 'Good triangle: Run %d - Run %d - Run %d'%(i, j, k)
    
    runOrientations = transformations_RR[0, :]
    outputFolder = './Output_runMerging'
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    joblib.dump(runOrientations, '%s/runOrientations.jbl'%outputFolder)
    

    # APPLY RUN TRANSFORMATIONS
    for i in range(0, nRuns):
        runNumber = runNumbers[i]
        runOrientation = runOrientations[i]
        print 'Run %s Orientation (wrt 0195) %s'%(runNumber, runOrientation)
        if not numpy.isnan(runOrientation):
            runOrientationMatrix = indexToMatrix(runOrientation)            
            spotMatricesList = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(runNumber, runNumber, runNumber)) 
            transformedSpotMatricesList = []
            
            for spotMatrix in spotMatricesList: # LATTICE:   n h k qRod I Icorrected h_transformed k_transformed Iscaled
                print 'spotsMatrix: (%d, %d)'%(spotMatrix.shape[0], spotMatrix.shape[1])
                if spotMatrix.shape[1] == 9:
                    transformedSpotsMatrix = []                    
                    for spot in spotMatrix:
                        h = spot[6]
                        k = spot[7]
                        indices = numpy.matrix('%d; %d'%(h, k)) #(2x1)
                        transformedIndices = runOrientationMatrix*indices
                        h_t = transformedIndices[0, 0]
                        k_t = transformedIndices[1, 0]
                        transformedSpot = [spot[0], h_t, k_t, spot[3], spot[8]]
                        transformedSpotsMatrix.append(transformedSpot)
                    transformedSpotsMatrix = numpy.asarray(transformedSpotsMatrix) # n h k qRod Iscaled
                    transformedSpotMatricesList.append(transformedSpotsMatrix)
                    print 'transformedSpotsMatrix: (%d, %d)'%(transformedSpotsMatrix.shape[0], transformedSpotsMatrix.shape[1])
        
        outputPath = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, runNumber)
        if not os.path.exists(outputPath):
            os.mkdir(outputPath)
        joblib.dump(transformedSpotMatricesList, '%s/r%s_transformedSpotsMatricesList.jbl'%(outputPath, runNumber))    
        
    
if __name__ == "__main__":
    print "\n**** CALLING transform_CCmethod_mergeRuns ****"
    main(sys.argv[1:])