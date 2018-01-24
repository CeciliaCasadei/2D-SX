# -*- coding: utf-8 -*-
import getopt
import joblib
import numpy
import sys
import os

import scaling_total_I


def scaling_mergeRunsFunction(myArguments):
    
    runNumbers = ['0195', '0196', '0197', '0200', '0201']
    
    # READ INPUTS  
    str_input = '--resolution_3D <resolution_3D> --cellSize <cellSize>'
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["resolution_3D=",
                                               "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python scaling_mergeRuns_total_I.py %s'%str_input
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_mergeRuns_total_I.py %s'%str_input
            sys.exit()
        elif option == "--resolution_3D":
            resolution_3D = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
    
    I_avg_singleRuns = []
    for runNumber in runNumbers:            
        inputPath = ('./Output_runMerging/spotsMatricesList-Transformed-r%s'
                     %(runNumber))        
        lattices = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                                %(inputPath, runNumber))
                                
        avg_Is = scaling_total_I.calculate_latticeAvgI(lattices, 
                                                       cellSize, 
                                                       resolution_3D)    
        avg_Is_cleaned = [avg_Is[j] for j in range(len(avg_Is)) 
                                    if not numpy.isnan(avg_Is[j])]   
        I_avg_singleRun = numpy.average(avg_Is_cleaned)  
        I_avg_singleRuns.append(I_avg_singleRun)
        print 'RUN: %s \nI AVG ON LATTICES: %.2f ph (AFTER SCALING)'%(runNumber,
                                                                      I_avg_singleRun)
    print 'SINGLE RUN AVERAGES:', I_avg_singleRuns
    tot_I_avg = numpy.average(I_avg_singleRuns)
    print 'TOTAL AVERAGE: ', tot_I_avg
    
    K_run = []
    for I_avg_singleRun in I_avg_singleRuns:
        K = tot_I_avg/I_avg_singleRun
        K_run.append(K)
    print 'SINGLE RUN SCALES: ', K_run                        

    # APPLY RUN SCALE FACTORS
    for runIndex in range(0, len(runNumbers)):
        runNumber = runNumbers[runIndex]
        runScale = K_run[runIndex] 
        inputPath = ('./Output_runMerging/spotsMatricesList-Transformed-r%s'
                     %(runNumber))        
        lattices = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                                %(inputPath, runNumber))
        scaledRun = []
        n = 0
        for lattice in lattices:
            lattice = numpy.asarray(lattice)
            scaledLattice = []
            if lattice[0, 4] == 1:  
                n = n+1
                for spot in lattice:
                    I_scaled = runScale * spot[3]
                    scaled_spot = [spot[0],  # h
                                   spot[1],  # k
                                   spot[2],  # qRod
                                   I_scaled, # I
                                   spot[4],  # flag
                                   spot[5],  # i_unassembled
                                   spot[6]]  # j_unassembled
                    scaledLattice.append(scaled_spot)
            else: 
                for spot in lattice:
                    scaled_spot = [spot[0],  # h
                                   spot[1],  # k
                                   spot[2],  # qRod
                                   spot[3],  # I
                                   spot[4],  # flag = 0
                                   spot[5],  # i_unassembled
                                   spot[6]]  # j_unassembled
                    scaledLattice.append(scaled_spot)
            scaledLattice = numpy.asarray(scaledLattice)
            scaledRun.append(scaledLattice)
        scaleOutputFolder = ('./Output_runMerging/spotsMatricesList-Scaled-r%s'
                             %runNumber)
        if not os.path.exists(scaleOutputFolder):
            os.mkdir(scaleOutputFolder)
        joblib.dump(scaledRun, 
                    '%s/r%s_scaledSpotsMatricesList.jbl'%(scaleOutputFolder, 
                                                          runNumber))       
        print runNumber, n, len(scaledRun)
        
           
if __name__ == "__main__":
    print "\n**** CALLING scaling_mergeRuns ****"
    scaling_mergeRunsFunction(sys.argv[1:])    