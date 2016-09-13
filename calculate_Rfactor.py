# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy

def calculate_Rfactor_Function(myArguments):
    # DEFAULTS
    inputFolder = ''
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    rodIndices = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                   
    
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

            
    sum_1 = 0
    sum_2 = 0
    
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector = []
        Irod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I) POINTS FROM ALL RUNS (AFTER RUN SCALING, WITH AVG SCALE SET TO 1)              
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(inputFolder, runNumber, runNumber))
        
            for latticeMatrix in myList:   # h k qRod I flag i_unassembled j_unassembled
                latticeMatrix = numpy.asarray(latticeMatrix)
                if latticeMatrix[0, 4] == 1:
                    for spot in latticeMatrix:
                        h = spot[0]
                        k = spot[1]
                        if (h == hRod and k == kRod) or (h == -hRod-kRod and k == hRod) or (h == kRod and k == -hRod-kRod):
                            Irod_vector.append(spot[3])
                            Qrod_vector.append(spot[2])
                        if (h == -hRod and k == -kRod) or (h == hRod+kRod and k == -hRod) or (h == -kRod and k == hRod+kRod):
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
                    sum_2 = sum_2 + I_mean
                if N > 1:
                    sum_reflection = sum_reflection * numpy.sqrt(float(N)/(N-1))
                else:
                    print 'Warning N = 1'
                    
                sum_1 = sum_1 + sum_reflection
                
    R_meas = float(sum_1)/sum_2
    print '\n\n'
    print R_meas
    fOpen = open('%s/R_meas.txt'%inputFolder, 'w')
    fOpen.write('R_meas = %.4f'%R_meas)
    fOpen.close()
                
        

if __name__ == "__main__":
    print "\n**** CALLING calculate_Rfactor ****"
    calculate_Rfactor_Function(sys.argv[1:])    