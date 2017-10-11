# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_merging.py
    
    
# DETERMINE AND APPLY RUN - RUN TRANSFORMATIONS
deltaQrodThreshold = 0.005    # A-1
n_minThreshold = 6
nLatticePairs = 1000

flag = 0
if flag == 1:
    os.system('python transform_CCmethod_mergeRuns.py --dQrod %f --nMin %d --nLatticePairs %d'
              %(deltaQrodThreshold, n_minThreshold, nLatticePairs))    
    


# DETERMINE AND APPLY RUN - RUN SCALE FACTORS
deltaQrodThreshold = 0.003    # A-1
n_minThreshold = 8
nLatticePairs = 1400
resolution_3D = 6.5           # A

flag = 0
if flag == 1:
    os.system('python scaling_mergeRuns.py --dQrod %f --nMin %d --nLatticePairs %d --resolution_3D %f'
              %(deltaQrodThreshold, n_minThreshold, nLatticePairs, resolution_3D))
    
    
# PLOT MERGED RODS AND DO POLYNOMIAL FIT
inputFolder = './Output_runMerging'
resolutionLimit = 7.0

flag = 0
if flag == 1:
    os.system('python merging.py --inputFolder %s --resolutionLimit %f'%(inputFolder, resolutionLimit))


newFolder = './Output_runMergingVsModel'
if not os.path.exists(newFolder):
    os.mkdir(newFolder)   