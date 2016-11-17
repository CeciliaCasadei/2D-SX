# -*- coding: utf-8 -*-
import os


# NORMALIZE LATTICES TO MODEL SCALE FACTOR AND APPLY SCALE FACTOR
flag = 0
if flag == 1:
    os.system('python model_applyScales.py')
    
    
# PLOT MERGED RODS
inputFolder = './Output_runMergingVsModel'
flag = 0
if flag == 1:
    os.system('python merging.py --inputFolder %s'%inputFolder)
    
    
# CALCULATE R-FACTOR
flag = 0
if flag == 1:
    os.system('python calculate_Rfactor.py --inputFolder %s'%inputFolder)


# CALCULATE PATTERSON
flag = 0
if flag == 1:
    os.system('python calculate_Patterson_normalized.py --inputFolder %s'%inputFolder)
    


# EXPORT h k l F sig(F) to MR
flag = 1
if flag == 1:
    os.system('python prepare_MRdata.py')
