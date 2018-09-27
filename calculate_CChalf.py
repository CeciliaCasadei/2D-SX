# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import random
import os

import simulate_resolution
import shannonSamplings
import get_rodIndices
import correlate

from bins import binLimits

# CALCULATE CChalf OF A RESOLUTION BIN [high, low] ONCE
def calculate(low, high, rodIndices, inputFolder, cellSize, overSampling, d):
    
    delta_qRod = (2*numpy.pi)/(overSampling*2*d) 
    
    N_uniques = 0
    CChalf_value = 0
    Is_1 = []
    Is_2 = []
    
    # Loop on Bragg lines
    for rod_hk in rodIndices:
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                      %(inputFolder, h, k))
        
        experimental_q = braggRodObject.experimental_q
        experimental_I = braggRodObject.experimental_I
        qMin = braggRodObject.qMin
        qMax = braggRodObject.qMax
         
        # Generate qRod Shannon values
        samplings = shannonSamplings.get_shannonSamplings(delta_qRod, qMax)
        samplings = samplings[1:-1]
        if not (samplings[0] > qMin and samplings[-1] < qMax):
            raise NameError('Samplings exceeding experimental range.')
            
        # Loop on qRod Shannon values
        for qRod in samplings:
            resolution = simulate_resolution.resolution(cellSize, h, k, qRod)
            if (low > resolution >= high):
                
                # Unique reflection belongs to resolution range
                N_uniques = N_uniques + 1
                
                shannonBin_l = qRod - (delta_qRod/2)
                shannonBin_r = qRod + (delta_qRod/2)

                Is_shannonBin = [experimental_I[i] 
                                 for i in range(0, len(experimental_q)) 
                                 if shannonBin_l < experimental_q[i] < shannonBin_r]
                N_obs = len(Is_shannonBin)
                
                Is_shannonBin_1 = []
                Is_shannonBin_2 = []
                
                # Extract half of the observations      
                random_sample = random.sample(range(N_obs), N_obs/2)
                for i in range(0, N_obs):
                    I_obs = Is_shannonBin[i]
                    if i in random_sample:
                        Is_shannonBin_1.append(I_obs)
                    else:
                        Is_shannonBin_2.append(I_obs)
                if not (N_obs == len(Is_shannonBin_1) + len(Is_shannonBin_2)):
                    raise Exception('Bad sub-group division.')
                    
                I1 = numpy.average(Is_shannonBin_1)
                I2 = numpy.average(Is_shannonBin_2)
                Is_1.append(I1)
                Is_2.append(I2)
                
                print ('SPOT %4d %4d %6.2f (%6.2f A) IN RANGE: %6.2f A - %6.2f A. AVG I1 AND I2: %6.2f %6.2f'
                       %(h, k, qRod, resolution, low, high, I1, I2))
                                
    CChalf_value = correlate.Correlate(Is_1, Is_2)
    
    if not (N_uniques == len(Is_1)):
        raise Exception('Error in N of unique reflections.')
        
    return N_uniques, CChalf_value

# CALCULATE CChalf OF A RESOLUTION BIN [high, low] n TIMES
def CChalf_f(low, 
           high, 
           rodIndices, 
           inputFolder, 
           cellSize, 
           overSampling, 
           d):
               
    n = 10
    CChalf_list = []
    N_uniques_list = []
    
    print '****** ', low, high, ' ******'
    
    for i in range(0, n):
        
        print '\n'
        N_uniques, CChalf_value = calculate(low, 
                                            high, 
                                            rodIndices, 
                                            inputFolder, 
                                            cellSize, 
                                            overSampling, 
                                            d)
        CChalf_list.append(CChalf_value)
        
        if len(N_uniques_list) > 0:
            if N_uniques != N_uniques_list[-1]:
                raise Exception('Error: different N_unique values.')
                
        N_uniques_list.append(N_uniques)
        
    return N_uniques_list, CChalf_list
    
def calculate_CCstar(CChalf):
    return numpy.sqrt( (2*CChalf)/(1+CChalf) )      

def calculate_CChalf_Function(myArguments):
    
    # SETTINGS
    nBins = len(binLimits) - 1
    resolution_2D = 6.0
    
    # KEYWORD ARGUMENTS
    input_str_0 = '--inputFolder <inputFolder>'
    input_str_1 = '--thickness <thickness> --overSampling <overSampling>'
    input_str_2 = '--cellSize <cellSize>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder=", 
                                                                 "thickness=",
                                                                 "overSampling=",
                                                                 "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python calculate_CChalf.py %s %s %s %s'%(input_str_0,
                                                               input_str_1,
                                                               input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_CChalf.py %s %s %s %s'%(input_str_0,
                                                                   input_str_1,
                                                                   input_str_2)
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
        elif option == "--thickness":
            d = float(value)
        elif option == "--overSampling":
            overSampling = int(value)
        elif option == "--cellSize":
            cellSize = float(value)
                
    # DEFINE ROD INDICES  
    rodIndices = get_rodIndices.defineRodIndices(resolution_2D)
    print '%d Rods to %.1f 2D-resolution'%(len(rodIndices), resolution_2D)         
    
    # LOG
    fOpen = open('%s/CChalf_bins_OS_%d.txt'%(inputFolder,
                                             overSampling), 'w')
    fOpen.write('Resolution 3D           N_uniques      CChalf   CCstar\n')
    
    # BINARY TO SAVE
    data = []
    
    # CALCULATE CChalf IN EACH 3D-RESOLUTION SHELL
    for i in range(0, nBins):
        low  = binLimits[i]
        high = binLimits[i+1]
        N_vector, CChalf_vector = CChalf_f(low, 
                                           high, 
                                           rodIndices, 
                                           inputFolder, 
                                           cellSize,
                                           overSampling, 
                                           d)
        CChalf_value = numpy.average(CChalf_vector)
        CCstar = calculate_CCstar(CChalf_value)
        N_uniques = N_vector[0]
        print CChalf_vector, CChalf_value, CCstar
        print N_vector
        fOpen.write('%6.2f - %6.2f      %12d       %.4f     %.4f\n'%(low, 
                                                                     high, 
                                                                     N_uniques, 
                                                                     CChalf_value,
                                                                     CCstar))
                                                             
        data_line = [low, high, N_uniques, CChalf_value, CCstar]
        data.append(data_line)
     
    # CALCULATE GLOBAL CChalf AND CCstar                                                        
    N_uniques_tot = 0
    for dl in data:
        N_uniques_tot = N_uniques_tot + dl[2]
    print "N_uniques_tot: ", N_uniques_tot
    
    avg_CChalf = 0
    avg_CCstar = 0
    for dl in data:
        weight = float(dl[2])/N_uniques_tot
        CChalf = dl[3]
        CCstar = dl[4]
        avg_CChalf = avg_CChalf + weight * CChalf
        avg_CCstar = avg_CCstar + weight * CCstar

    fOpen.write('%6.2f - %6.2f      %12d       %.4f     %.4f\n'%(binLimits[0], 
                                                                 binLimits[-1], 
                                                                 N_uniques_tot, 
                                                                 avg_CChalf,
                                                                 avg_CCstar))
    
    data_line = [binLimits[0], binLimits[-1], N_uniques_tot, avg_CChalf, avg_CCstar]
    data.append(data_line)
                                                                      
    fOpen.close()
    
    data = numpy.asarray(data)
    dataBinary_folder = '%s/CChalf'%inputFolder
    if not os.path.exists(dataBinary_folder):
        os.mkdir(dataBinary_folder)
    joblib.dump(data, '%s/CChalf_bins_OS_%d.jbl'%(dataBinary_folder,
                                                  overSampling))

if __name__ == "__main__":
    print "\n**** CALLING calculate_CChalf ****"
    calculate_CChalf_Function(sys.argv[1:])    
    