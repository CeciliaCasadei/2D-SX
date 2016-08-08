# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy

def calculate_Rfactor_Function(myArguments):
    
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h")
    except getopt.GetoptError:
        print 'Usage: python plotRods.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRods.py'
            sys.exit()

        
    rodIndices = [[1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                   
    
    sum_1 = 0
    sum_2 = 0
    # PLOT I vs QROD, FOR EVERY ROD, AND PRODUCE MODEL (POLYNOMIAL FIT OF MEDIAN VALUES)    
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector = []
        Irod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I) POINTS FROM ALL RUNS (AFTER RUN SCALING, WITH AVG SCALE SET TO 1)              
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            myList = joblib.load('./Output_runMerging/spotsMatricesList-Scaled-r%s-AvgTo1/r%s_scaledSpotsMatricesList.jbl'%(runNumber, runNumber))
        
            for latticeMatrix in myList:   # n h k qRod Iscaled_old I_scaled     
                for spot in latticeMatrix:
                    h_transformed = spot[1]
                    k_transformed = spot[2]
                    if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                        Irod_vector.append(spot[5])
                        Qrod_vector.append(spot[3])
                    if (h_transformed == -hRod and k_transformed == -kRod) or (h_transformed == hRod+kRod and k_transformed == -hRod) or (h_transformed == -kRod and k_transformed == hRod+kRod):
                        Irod_vector.append(spot[5])
                        Qrod_vector.append(-spot[3])
         
        # REMOVE NAN VALUES
        cleanedList_Irod = [Irod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
        cleanedList_Qrod = [Qrod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
                    
        # BINNING
        Qrod_min = min(cleanedList_Qrod)
        Qrod_max = max(cleanedList_Qrod)
        bins, step = numpy.linspace(Qrod_min, Qrod_max, num=(Qrod_max-Qrod_min)/0.002, endpoint = True, retstep = True)
        
        
        # LOOP ON BINS
        for i in range(0, len(bins)-1):
            edge_l = bins[i]
            edge_r = bins[i+1]
            
            binList_Irod = [cleanedList_Irod[item] for item in range(0, len(cleanedList_Irod)) if edge_l <= cleanedList_Qrod[item] <= edge_r]
            N = len(binList_Irod)
            if N > 0:
                
                I_mean = numpy.average(binList_Irod)
                sum_reflection = 0
                for I in binList_Irod:
                    sum_reflection = sum_reflection + abs(I - I_mean)
                    sum_2 = sum_2 + I
                if N > 1:
                    sum_reflection = sum_reflection * numpy.sqrt(float(N)/(N-1))
                else:
                    print 'Warning N = 1'
                    
                sum_1 = sum_1 + sum_reflection
                
        R_meas = float(sum_1)/sum_2
    print '**********************************************************************************************'
    print R_meas
                
                
                
def calculate_Rfactor_Function_model(myArguments):
    
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h")
    except getopt.GetoptError:
        print 'Usage: python plotRods.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRods.py'
            sys.exit()
        


    
    rodIndices = [[1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                  
    sum_1 = 0
    sum_2 = 0

    # PLOT I vs QROD, FOR EVERY ROD, AND PRODUCE MODEL (POLYNOMIAL FIT OF MEDIAN VALUES)    
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector = []
        Irod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I) POINTS FROM ALL RUNS            
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            myList = joblib.load('./Output_runMerging/transformAndScaleToModel_r%s/r%s_scaledLatticesList/r%s_scaledLatticesList.jbl'%(runNumber, runNumber, runNumber))
            #myList = joblib.load('./Output_runMerging/spotsMatricesList-Scaled-r%s-AvgTo1/r%s_scaledSpotsMatricesList.jbl'%(runNumber, runNumber))
        
            for latticeMatrix in myList:   # h_transformed k_transformed qRod I_scaled   
                for spot in latticeMatrix:
                    h_transformed = spot[0]
                    k_transformed = spot[1]
                    if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                        Irod_vector.append(spot[3])
                        Qrod_vector.append(spot[2])
                    if (h_transformed == -hRod and k_transformed == -kRod) or (h_transformed == hRod+kRod and k_transformed == -hRod) or (h_transformed == -kRod and k_transformed == hRod+kRod):
                        Irod_vector.append(spot[3])
                        Qrod_vector.append(-spot[2])
         
        # REMOVE NAN VALUES
        cleanedList_Irod = [Irod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
        cleanedList_Qrod = [Qrod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
                    
        # BINNING
        Qrod_min = min(cleanedList_Qrod)
        Qrod_max = max(cleanedList_Qrod)
        bins, step = numpy.linspace(Qrod_min, Qrod_max, num=(Qrod_max-Qrod_min)/0.002, endpoint = True, retstep = True)
        
        # LOOP ON BINS
        for i in range(0, len(bins)-1):
            edge_l = bins[i]
            edge_r = bins[i+1]
            
            binList_Irod = [cleanedList_Irod[item] for item in range(0, len(cleanedList_Irod)) if edge_l <= cleanedList_Qrod[item] <= edge_r]
            N = len(binList_Irod)
            if N > 0:
                
                I_mean = numpy.average(binList_Irod)
                sum_reflection = 0
                for I in binList_Irod:
                    sum_reflection = sum_reflection + abs(I - I_mean)
                    sum_2 = sum_2 + I
                if N > 1:
                    sum_reflection = sum_reflection * numpy.sqrt(float(N)/(N-1))
                else:
                    print 'Warning N = 1'
                    
                sum_1 = sum_1 + sum_reflection
                
    R_meas = float(sum_1)/sum_2
    print R_meas
                
        
                
        
        

if __name__ == "__main__":
    print "\n**** CALLING calculate_Rfactor ****"
    calculate_Rfactor_Function(sys.argv[1:])    
    calculate_Rfactor_Function_model(sys.argv[1:])    
    