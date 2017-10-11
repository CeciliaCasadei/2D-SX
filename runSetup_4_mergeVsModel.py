# -*- coding: utf-8 -*-
import os


# NORMALIZE LATTICES TO MODEL SCALE FACTOR AND APPLY SCALE FACTOR
flag = 0
if flag == 1:
    os.system('python model_applyScales.py')
    
    
# PLOT MERGED RODS, DO POLYNOMIAL FIT.
inputFolder = './Output_runMergingVsModel'
resolutionLimit = 6.0

flag = 0
if flag == 1:
    os.system('python merging.py --inputFolder %s --resolutionLimit %f'%(inputFolder, resolutionLimit))