# -*- coding: utf-8 -*-
import joblib
import os
import numpy
import sys
import getopt

import makeOrbits

def sinc_function(x, a, N, delta, T, d):
    y = 0
    for i in range(-N, N+1):
        j = i+N
        if abs(x-(i*delta)) < 0.001: # Possible only if x is a single value
            y = y + a[j]
        else:
            y = (y + 
                 a[j]  * (numpy.sin(d*(x-i*delta)) / (d*(x-i*delta)) ) 
                       * numpy.exp(-0.5*(1/(4*numpy.pi**2))*T*(x-i*delta)**2))
    return y
    
def prepare_MRdata_Function(myArguments):
    input_str_1 = '--resolutionLimit <resolutionLimit>'
    input_str_2 = '--thickness <thickness> --damping <damping>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["resolutionLimit=",
                                               "thickness=",
                                               "damping="])
    except getopt.GetoptError:
        print 'Usage: python prepare_MRdata.py %s'%(input_str_1, 
                                                    input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python prepare_MRdata.py %s'%(input_str_1, 
                                                        input_str_2)
            sys.exit()
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--thickness":
            d = float(value)
        elif option == "--damping":
            T = float(value)

    outFolder = './Output_runMergingVsModel/dataToMR'
    if not os.path.exists(outFolder):
        os.mkdir(outFolder)
    
    hkl_file = open('%s/experimentalReflections_h_k_l_F_sigF.hkl'
                     %outFolder,'w')
    hkl_file_q_I = open('%s/experimentalReflections_h_k_l_qRod_I_sigI.hkl'
                         %outFolder, 'w')
    
    d = 45
    c_star = (2*numpy.pi)/(2*d)
    T = 80
       
    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices_all = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices_all.append(orbit_label)        
    print '%d Rods'%len(rodIndices_all) 

    # LOOP ON BRAGG RODS
    for rodIndices in rodIndices_all:
        print rodIndices
        h = rodIndices[0]
        k = rodIndices[1]
        
        braggRodObjectFile = ('./Output_runMergingVsModel/Shannon_sampling/braggRodObjects/braggRodObject_%d_%d.jbl'
                              %(h, k))
        braggRodObject = joblib.load(braggRodObjectFile)
        
        qRod_max = braggRodObject.qMax
        qRod_min = braggRodObject.qMin
        
        experimental_qs = braggRodObject.experimental_q
        experimental_Is = braggRodObject.experimental_I
        
        model_coefficients = braggRodObject.model_coefficients
        
        
        
        l_max = min( [int(qRod_max/c_star), int(-qRod_min/c_star)] )
        
    
        # DIVIDE EACH ROD IN INTERVALS CENTERED IN l VALUES
        for l in range(-l_max, l_max+1):
            
            
            qRod_center = l*c_star
            print 'qRod_center ', qRod_center
            
            qRod_left   = qRod_center - (c_star/2)
            qRod_right  = qRod_center + (c_star/2)
            
            if qRod_center < qRod_min or qRod_center > qRod_max:
                print 'PROBLEM l=%d'%l
            
            I_model_center = sinc_function(qRod_center, 
                                           model_coefficients, 
                                           (len(model_coefficients)-1)/2, 
                                           c_star, 
                                           T, 
                                           d)
            if I_model_center >= 0:
                F = numpy.sqrt(I_model_center)
                N = 0
                sumOfSquaredResiduals = 0
                for index in range(0, len(experimental_qs)):
                    if qRod_left < experimental_qs[index] < qRod_right:
                        N = N + 1
                        q_experimental = experimental_qs[index]
                        I_experimental = experimental_Is[index]
                        I_model = sinc_function(q_experimental, 
                                                model_coefficients, 
                                                (len(model_coefficients)-1)/2, 
                                                c_star, 
                                                T, 
                                                d)
                        squaredResidual = (I_experimental-I_model)**2
                        sumOfSquaredResiduals = sumOfSquaredResiduals + squaredResidual
                    
                varI = sumOfSquaredResiduals/N
                sigI = numpy.sqrt(varI)
                sigF = numpy.sqrt(I_model_center + sigI) - F
                hkl_file.write('%3d %3d %3d %8.2f %8.2f\n'%(h, k, l, F, sigF))
                hkl_file_q_I.write('%3d %3d %3d %8.2f %8.2f %8.2f\n'
                                    %(h, k, l, qRod_center, I_model_center, sigI))
            else:
                F = 0
                sigF = 0
            
            print '%3d %3d %3d %8.2f %8.2f'%(h, k, l, F, sigF)
            
    hkl_file.close()
    hkl_file_q_I.close()

if __name__ == "__main__":
    print "\n**** CALLING prepare_MRdata ****"
    prepare_MRdata_Function(sys.argv[1:])    
            