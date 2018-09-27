# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os

import shannon_model
import simulate_resolution
import get_rodIndices
import shannonSamplings
import FW_functions

from bins import binLimits



def calculate_binAvgI(low, 
                      high, 
                      rodIndices, 
                      inputFolder, 
                      cellSize, 
                      T, 
                      d, 
                      overSampling):
    
    c_star = (2*numpy.pi)/(2*d) 
    Is_bin = []
    
    for rod_hk in rodIndices:
        
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                      %(inputFolder, h, k))
    
        qMax = braggRodObject.qMax
        model_coefficients = braggRodObject.model_coefficients
        
        samplingSpace = c_star/overSampling  
        samplings = shannonSamplings.get_shannonSamplings(samplingSpace, qMax)
        samplings = samplings[1:-1]
        
        for q in samplings:
            resolution = simulate_resolution.resolution(cellSize, 
                                                        h, 
                                                        k, 
                                                        q)
            if high <= resolution < low:
                # VALUE OF SHANNON MODEL I IN THE MIDDLE OF SHANNON BIN
                I = shannon_model.sinc_function(q, 
                                                model_coefficients, 
                                                (len(model_coefficients)-1)/2, 
                                                c_star, 
                                                T, 
                                                d)
                Is_bin.append(I)
      
    Is_bin = numpy.asarray(Is_bin)
    binAvgI = numpy.average(Is_bin)    
    return binAvgI
    
    

def FW(low, high, rodIndices, inputFolder, cellSize, T, d, overSampling):
    
    I_FW_over_sigI_FW_vector = []
    data = []
    c_star = (2*numpy.pi)/(2*d) 
    
    binAvgI = calculate_binAvgI(low, 
                                high, 
                                rodIndices, 
                                inputFolder, 
                                cellSize, 
                                T, 
                                d, 
                                overSampling)
     
    print '*** ', low, high, ' *** ' , binAvgI, ' ***'
    for rod_hk in rodIndices:
        
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                      %(inputFolder, h, k))
        
        experimental_q = braggRodObject.experimental_q
        experimental_I = braggRodObject.experimental_I
        qMax = braggRodObject.qMax
        model_coefficients = braggRodObject.model_coefficients
               
        ##########################
        samplingSpace = c_star/overSampling  
        ##########################
        
        samplings = shannonSamplings.get_shannonSamplings(samplingSpace, qMax)
        samplings = samplings[1:-1]
        
        for q in samplings:
            resolution = simulate_resolution.resolution(cellSize, 
                                                        h, 
                                                        k, 
                                                        q)
            if high <= resolution < low:
                # VALUE OF SHANNON MODEL I IN THE MIDDLE OF SHANNON BIN
                I = shannon_model.sinc_function(q, 
                                                model_coefficients, 
                                                (len(model_coefficients)-1)/2, 
                                                c_star, 
                                                T, 
                                                d)
                                
                # STD DEVIATION WITHIN SHANNON BIN
                qLeft  = q - samplingSpace/2
                qRight = q + samplingSpace/2
                
                qs_bin = [experimental_q[i] for i in range(0, len(experimental_q))
                                            if qLeft < experimental_q[i] < qRight]
                Is_bin = [experimental_I[i] for i in range(0, len(experimental_q))
                                            if qLeft < experimental_q[i] < qRight]
                                                
                sum_delta_sq = 0
                N = 0
                for i in range(0, len(qs_bin)):
                    qobs = qs_bin[i]
                    Iobs = Is_bin[i]
                    Imodel = shannon_model.sinc_function(qobs, 
                                                         model_coefficients, 
                                                         (len(model_coefficients)-1)/2, 
                                                         c_star, 
                                                         T, 
                                                         d)
                    delta_sq = (Iobs-Imodel)**2
                    sum_delta_sq = sum_delta_sq + delta_sq
                    N = N+1
                    
                # Formula for the std dev of sample mean
                sigma_I = numpy.sqrt(sum_delta_sq)/N   
                
                correction = (sigma_I**2)/binAvgI
                I_corrected = I - correction
                
                if I>-4*sigma_I:         
                    I_FW    = FW_functions.calc_avg_intensity(I_corrected, 
                                                              sigma_I)
                    sigI_FW = numpy.sqrt(FW_functions.calc_var_intensity(I_corrected, 
                                                                         sigma_I))
                    F_FW    = FW_functions.calc_avg_amplitude(I_corrected, 
                                                              sigma_I)
                    sigF_FW = numpy.sqrt(FW_functions.calc_var_amplitude(I_corrected, 
                                                                         sigma_I))
                    I_FW_over_sigI_FW = I_FW / sigI_FW
                    I_FW_over_sigI_FW_vector.append(I_FW_over_sigI_FW)
                    data_line = [h, 
                                 k, 
                                 int(round(q/samplingSpace)),
                                 q, 
                                 I, 
                                 sigma_I, 
                                 I_FW, 
                                 sigI_FW, 
                                 F_FW, 
                                 sigF_FW,
                                 I_FW_over_sigI_FW]
                    data.append(data_line)
                else:
                    print 'NAN'
                    I_FW = numpy.nan
                    sigI_FW = numpy.nan
                    F_FW = numpy.nan
                    sigF_FW = numpy.nan
                    I_FW_over_sigI_FW = numpy.nan
                    data_line = [h, 
                                 k, 
                                 int(round(q/samplingSpace)),
                                 q, 
                                 I, 
                                 sigma_I, 
                                 I_FW, 
                                 sigI_FW, 
                                 F_FW, 
                                 sigF_FW,
                                 I_FW_over_sigI_FW]
                    data.append(data_line)
                                                                                                                                                
    I_FW_over_sigI_FW_vector = numpy.asarray(I_FW_over_sigI_FW_vector)
    StoN = numpy.average(I_FW_over_sigI_FW_vector)
    data = numpy.asarray(data)
    
    return len(I_FW_over_sigI_FW_vector), StoN, data

        

def calculate_FW_Function(myArguments):
    
    # SETTINGS
    nBins = len(binLimits) - 1
    resolution_2D = 6.0
    
    # KEYWORD ARGUMENTS
    input_str_0 = '--inputFolder <inputFolder>'
    input_str_1 = '--thickness <thickness> --damping <damping>'
    input_str_2 = '--cellSize <cellSize> --overSampling <overSampling>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder=", 
                                                                 "thickness=",
                                                                 "damping=",
                                                                 "cellSize=",
                                                                 "overSampling="])
    except getopt.GetoptError:
        print 'Usage: python calculate_FW_values.py %s %s %s %s'%(input_str_0,
                                                                  input_str_1,
                                                                  input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python calculate_FW_values.py %s %s %s %s'%(input_str_0,
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
        elif option == "--overSampling":
            overSampling = int(value)
    
    # SAVE BINARIES
    dataBinary_folder = '%s/French_Wilson_OS_%d'%(inputFolder, overSampling)
    if not os.path.exists(dataBinary_folder):
        os.mkdir(dataBinary_folder)        
        
    # DEFINE ROD INDICES  
    rodIndices = get_rodIndices.defineRodIndices(resolution_2D)
    print '%d Rods to %.1f 2D-resolution'%(len(rodIndices), resolution_2D) 
    
    # LOG
    fOpen = open('%s/I_over_sigI_FW_bins_OS_%d.txt'%(inputFolder,
                                                     overSampling), 'w')
    fOpen.write('Resolution 3D           N            I_FW / sigI_FW \n')
    
    # BINARY TO SAVE
    data_bins = []
    
    # CALCULATE FW VALUES IN EACH 3D-RESOLUTION SHELL
    for i in range(0, nBins):
        low  = binLimits[i]
        high = binLimits[i+1]
        N, StoN, data = FW(low, 
                           high, 
                           rodIndices, 
                           inputFolder, 
                           cellSize,
                           T, 
                           d,
                           overSampling)
        
        fOpen.write('%6.2f - %6.2f      %12d       %.2f \n'%(low, 
                                                             high, 
                                                             N, 
                                                             StoN))
        data_bin = [low, high, N, StoN]
        data_bins.append(data_bin)
        
        joblib.dump(data, '%s/FW_uniques_bin_%d.jbl'%(dataBinary_folder, i))
               
    fOpen.close()

    data_bins = numpy.asarray(data_bins)    
    joblib.dump(data_bins, '%s/FW_bins.jbl'%dataBinary_folder)
    
    
        
def calculate_global(inputFolder, overSampling):
    nBins = len(binLimits) - 1
    dataBinary_folder = '%s/French_Wilson_OS_%d'%(inputFolder, overSampling)
    
    StoN_list = []
    for i in range(0,nBins):
        print i
        data = joblib.load('%s/FW_uniques_bin_%d.jbl'%(dataBinary_folder, i))
        for dataLine in data:
            StoN = dataLine[10]
            StoN_list.append(StoN)
                
    StoN_list = numpy.asarray(StoN_list)
    globalAvgStoN = numpy.average(StoN_list)
    print "OVERSAMPLING: ", overSampling
    print " N UNIQUES: ", len(StoN_list)
    print "I_FW/sig(I_FW): ", globalAvgStoN
    
    
def write_unique(inputFolder, overSampling):
    nBins = len(binLimits) - 1
    dataBinary_folder = '%s/French_Wilson_OS_%d'%(inputFolder, overSampling)
    
    
    fOpenFile = '%s/FW_uniques_OS_%d.txt'%(dataBinary_folder, overSampling)
    fOpen = open(fOpenFile, 'w')
    fOpen.write('    h     k     l    qz    Iobs sigIobs    I_FW sigI_FW    F_FW sigF_FW I_FW_over_sigI_FW\n')
        
    for i in range(0,nBins):
        print i
        data = joblib.load('%s/FW_uniques_bin_%d.jbl'%(dataBinary_folder, i))
        for data_line in data:
            fOpen.write('%5d %5d %5d %5.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n'%(data_line[0],
                                                                                         data_line[1],
                                                                                         data_line[2],
                                                                                         data_line[3],
                                                                                         data_line[4],
                                                                                         data_line[5],
                                                                                         data_line[6],
                                                                                         data_line[7],
                                                                                         data_line[8],
                                                                                         data_line[9],
                                                                                         data_line[10]))                                                                                                                        
    fOpen.close()
    
    
    
if __name__ == "__main__":
    print "\n**** CALLING calculate_FW_values ****"
    calculate_FW_Function(sys.argv[1:])  