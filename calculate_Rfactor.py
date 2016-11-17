# -*- coding: utf-8 -*-
import sys
import getopt
import joblib

def calculate_Rfactor_Function(myArguments):
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder="])
    except getopt.GetoptError:
        print 'Usage: python calculate_Rfactor.py --inputFolder <inputFolder>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_Rfactor.py --inputFolder <inputFolder>'
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value

    rodIndices = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]  
    
    R_value_up = 0
    R_value_down = 0    
    for rod_hk in rodIndices:
        
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'%(inputFolder, h, k))
        
        experimental_q = braggRodObject.experimental_q
        experimental_I = braggRodObject.experimental_I
        model_coefficients = braggRodObject.model_coefficients
        
        for index in range(0, len(experimental_q)):
            q_observed = experimental_q[index]
            I_observed = experimental_I[index]
            
            I_model = 0
            for exponent in range(0, len(model_coefficients)):
                I_model = I_model + model_coefficients[exponent] * (q_observed**exponent)
            R_value_up = R_value_up + abs(I_observed-I_model)
            R_value_down = R_value_down + abs(I_model)
            
    R_value = R_value_up/R_value_down
    print R_value
    fOpen = open('%s/R_factor.txt'%inputFolder, 'w')
    fOpen.write('R_factor = %.4f'%R_value)
    fOpen.close()

if __name__ == "__main__":
    print "\n**** CALLING calculate_Rfactor ****"
    calculate_Rfactor_Function(sys.argv[1:])    
    