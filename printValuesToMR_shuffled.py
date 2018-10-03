# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import random
from bins import binLimits



def printValues(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["overSampling="])
    except getopt.GetoptError:
        print 'Usage: python printValuesToMR.py --overSampling <OS>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python printValuesToMR.py --overSampling <OS>'
            sys.exit()
        elif option == "--overSampling":
            overSampling = value
    
    folder = './Output_runMergingVsModel/Shannon_sampling'
    
    nBins = len(binLimits) - 1
   
    fOpen = open('%s/h_k_l_F_sigF_FW_OS_%s_shuffled.txt'%(folder,
                                                          overSampling), 'w')
    
    for i in range(0, nBins):
        FW_data = joblib.load('%s/French_Wilson_OS_%s/FW_uniques_bin_%d.jbl'%(folder, 
                                                                              overSampling, 
                                                                              i))
        N = FW_data.shape[0]
        print FW_data.shape   
           
        for spot in FW_data:
            
            h    = int(spot[0])
            k    = int(spot[1])
            l    = int(spot[2])
            n    = random.randint(0, N-1)
            F    = FW_data[n, 8]
            sigF = FW_data[n, 9]
            
            fOpen.write('%5d %5d %5d %8.4f %8.4f\n'%(h,
                                                     k,
                                                     l,
                                                     F,
                                                     sigF))
           
    fOpen.close()

if __name__ == "__main__":
    print "\n**** CALLING printValuesToMR ****"
    printValues(sys.argv[1:])  