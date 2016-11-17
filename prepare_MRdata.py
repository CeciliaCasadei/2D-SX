# -*- coding: utf-8 -*-
import joblib
import os
import numpy
import sys
import getopt

def prepare_MRdata_Function(myArguments):
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h")
    except getopt.GetoptError:
        print 'Usage: python prepare_MRdata.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python prepare_MRdata.py'
            sys.exit()

    if not os.path.exists('./Output_data_to_MR'):
        os.mkdir('./Output_data_to_MR')
    
    hkl_file = open('./Output_data_to_MR/experimentalReflections_c_504p5A_F_sigF.hkl', 'w')
    hkl_file_q_I = open('./Output_data_to_MR/experimentalReflections_c_504p5A_qRod_I.hkl', 'w')
    
    c = 100.9
    c_long = 5*c
    c_long_star = 2*numpy.pi/c_long
    
    
    rodIndices_all = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                      [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                      [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                      
    # LOOP ON BRAGG RODS
    for rodIndices in rodIndices_all:
        print rodIndices
        h = rodIndices[0]
        k = rodIndices[1]
        
        braggRodObject = joblib.load('Output_runMergingVsModel/braggRodObjects/braggRodObject_%d_%d.jbl'%(h, k))
        
        experimental_qs = braggRodObject.experimental_q
        experimental_Is = braggRodObject.experimental_I
        model_coefficients = braggRodObject.model_coefficients
    
        qRod_max = braggRodObject.qMax
        qRod_min = braggRodObject.qMin
        
        l_max = min( [int(qRod_max/c_long_star), int(-qRod_min/c_long_star)] )
    
        # DIVIDE EACH ROD IN INTERVALS CENTERED IN l VALUES
        for l in range(-l_max, l_max+1):
            N = 0
            sumOfSquaredResiduals = 0
            
            qRod_center = l*c_long_star        
            qRod_left  = qRod_center - (c_long_star/2)
            qRod_right = qRod_center + (c_long_star/2)
            
            if qRod_center < qRod_min or qRod_center > qRod_max:
                print 'PROBLEM l=%d'%l
            
            I_model_center = 0
            for exponent in range(0, len(model_coefficients)):
                I_model_center = I_model_center + model_coefficients[exponent] * (qRod_center**exponent)
            if I_model_center > 0:
                F = numpy.sqrt(I_model_center)
            else:
                continue
            
            
                
            for index in range(0, len(experimental_qs)):
                if qRod_left < experimental_qs[index] < qRod_right:
                    N = N + 1
                    q_experimental = experimental_qs[index]
                    I_experimental = experimental_Is[index]
                    I_model = 0
                    for exponent in range(0, len(model_coefficients)):
                        I_model = I_model + model_coefficients[exponent] * (q_experimental**exponent)
                    squaredResidual = (I_experimental-I_model)**2
                    sumOfSquaredResiduals = sumOfSquaredResiduals + squaredResidual
                    
            varI = sumOfSquaredResiduals/N
            sigI = numpy.sqrt(varI)
            sigF = numpy.sqrt(I_model_center + sigI) - F
            hkl_file.write('%3d %3d %3d %8.2f %8.2f\n'%(h, k, l, F, sigF))
            hkl_file_q_I.write('%3d %3d %3d %8.2f %8.2f\n'%(h, k, l, qRod_center, I_model_center))
            print '%3d %3d %3d %8.2f %8.2f'%(h, k, l, F, sigF)
            
    hkl_file.close()

# 1466 reflections out of 4210    

if __name__ == "__main__":
    print "\n**** CALLING prepare_MRdata ****"
    prepare_MRdata_Function(sys.argv[1:])    
            