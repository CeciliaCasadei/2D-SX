# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_merging.py
    
    
# DETERMINE AND APPLY RUN - RUN TRANSFORMATIONS
deltaQrodThreshold = 0.005
n_minThreshold = 6
nLatticePairs = 200

flag = 0
if flag == 1:
    os.system('python transform_CCmethod_mergeRuns.py --dQrod %f --nMin %d --nLatticePairs %d'%(deltaQrodThreshold, n_minThreshold, nLatticePairs))    
    #os.system('python transform_CCmethod_mergeRuns.py --dQrod 0.005 --nMin 6 --nLatticePairs 200')
    


# DETERMINE AND APPLY RUN - RUN SCALE FACTORS
deltaQrodThreshold = 0.003
n_minThreshold = 8
nLatticePairs = 1400

flag = 0
if flag == 1:
    os.system('python scaling_mergeRuns.py --dQrod %f --nMin %d --nLatticePairs %d'%(deltaQrodThreshold, n_minThreshold, nLatticePairs))
    
    
# PLOT MERGED RODS
inputFolder = './Output_runMerging'
flag = 1
if flag == 1:
    os.system('python merging.py --inputFolder %s'%inputFolder)
    
    
# CALCULATE R-MEAS
flag = 1
if flag == 1:
    os.system('python calculate_Rfactor.py --inputFolder %s'%inputFolder)