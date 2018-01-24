# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy

import shannon_model
import simulate_resolution
import get_rodIndices


def R(low, high, rodIndices, inputFolder, cellSize, T, d):
    
    c_star = (2*numpy.pi)/(2*d) 
    
    R_value_up = 0
    R_value_down = 0   
    N = 0
    
    print '****** ', low, high, ' ******'
    for rod_hk in rodIndices:
        print rod_hk
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                      %(inputFolder, h, k))
        
        experimental_q = braggRodObject.experimental_q
        experimental_I = braggRodObject.experimental_I
        model_coefficients = braggRodObject.model_coefficients
        
        for index in range(0, len(experimental_q)):
            q_observed = experimental_q[index]
            
            resolution = simulate_resolution.resolution(cellSize, 
                                                        h, 
                                                        k, 
                                                        q_observed)
            if high <= resolution < low:
                I_observed = experimental_I[index]
                
                I_model = shannon_model.sinc_function(q_observed, 
                                                      model_coefficients, 
                                                      (len(model_coefficients)-1)/2, 
                                                      c_star, 
                                                      T, 
                                                      d)
                R_value_up = R_value_up + abs(I_observed-I_model)
                R_value_down = R_value_down + abs(I_model)
                N = N+1
            
    R_value = R_value_up/R_value_down
    return N, R_value
        

def calculate_Rfactor_Function(myArguments):
    
    # SETTINGS
    bins = [54.09, 25.27, 16.66, 12.79, 11.64, 10.19, 9.41, 8.64, 7.94, 7.53, 
            7.13, 6.79, 6.54, 6.20, 6.02, 5.30]
    nBins = len(bins) - 1
    resolution_2D = 6.0
    
    # KEYWORD ARGUMENTS
    input_str_0 = '--inputFolder <inputFolder>'
    input_str_1 = '--thickness <thickness> --damping <damping>'
    input_str_2 = '--cellSize <cellSize>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder=", 
                                                                 "thickness=",
                                                                 "damping=",
                                                                 "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python calculate_Rfactor.py %s %s %s %s'%(input_str_0,
                                                                input_str_1,
                                                                input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_Rfactor.py %s %s %s %s'%(input_str_0,
                                                                    input_str_1,
                                                                    input_str_2)
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
        elif option == "--thickness":
            d = float(value)
        elif option == "--damping":
            T = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
         
    # DEFINE ROD INDICES  
    rodIndices = get_rodIndices.defineRodIndices(resolution_2D)
    print '%d Rods to %.1f 2D-resolution'%(len(rodIndices), resolution_2D) 
    
    # LOG
    fOpen = open('%s/R_factor_bins.txt'%inputFolder, 'w')
    fOpen.write('Resolution 3D           N            R\n')
    
    # CALCULATE R IN EACH 3D-RESOLUTION SHELL
    for i in range(0, nBins):
        low = bins[i]
        high = bins[i+1]
        N, R_value = R(low, 
                       high, 
                       rodIndices, 
                       inputFolder, 
                       cellSize,
                       T, 
                       d)
        print R_value
        fOpen.write('%6.2f - %6.2f      %12d       %.2f \n'%(low, 
                                                             high, 
                                                             N, 
                                                             R_value))
                                                             
                                                             
    # CALCULATE GLOBAL R WITH DIFFERENT HIGH-RES CUTOFFS                                                         
    for secondEdge in [bins[-1], bins[-2]]:
        
        N, R_value = R(bins[0], 
                       secondEdge, 
                       rodIndices, 
                       inputFolder, 
                       cellSize, 
                       T, 
                       d)
        print 'Global: ', R_value
        fOpen.write('\n%6.2f - %6.2f      %12d       %.3f \n'%(bins[0], 
                                                               secondEdge, 
                                                               N, 
                                                               R_value))
                                                               
    
         
    fOpen.close()

if __name__ == "__main__":
    print "\n**** CALLING calculate_Rfactor ****"
    calculate_Rfactor_Function(sys.argv[1:])    
    