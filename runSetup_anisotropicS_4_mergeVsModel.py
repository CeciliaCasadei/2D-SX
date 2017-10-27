# -*- coding: utf-8 -*-
import os


# NORMALIZE LATTICES TO MODEL SCALE FACTOR AND APPLY SCALE FACTOR
cellSize = 62.45

flag = 0
if flag == 1:
    os.system('python model_applyScales_anisotropic.py --cellSize %f'%cellSize)
 
   
    
# PLOT MERGED RODS, DO POLYNOMIAL FIT.
inputFolder = './Output_runMergingVsModel'
resolutionLimit = 6.0

flag = 0
if flag == 1:
    os.system('python merging.py --inputFolder %s --resolutionLimit %f'%(inputFolder, resolutionLimit))
 
   

# PLOT MERGED RODS, BUILD SINC MODEL
resolutionLimit = 6.0      # A, 2D
thickness = 45             # A
damping = 80

flag = 1
if flag == 1:
    os.system('python rodsFit_shannonTheo.py \
               --resolutionLimit %f \
               --thickness %f \
               --damping %f \
               --folder %s'
               %(resolutionLimit, thickness, damping, inputFolder))
               


# COUNT N OBSERVATIONS VS Q_3D    
resolutionLimit = 6.0      # A, 2D

flag = 0
if flag == 1:
    os.system('python nObs.py \
               --resolutionLimit %f \
               --cellSize %f'
               %(resolutionLimit, cellSize))



# EXPORT h k l F sig(F) to MR
flag = 0
if flag == 1:
    os.system('python prepare_MRdata.py --resolutionLimit %f'%resolutionLimit)

    
    
## CALCULATE R-FACTOR
#flag = 0
#if flag == 1:
#    os.system('python calculate_Rfactor.py --inputFolder %s --resolutionLimit %f'%(inputFolder, resolutionLimit))
#
#
## CALCULATE PATTERSON
#flag = 0
#if flag == 1:
#    os.system('python calculate_Patterson_normalized.py --inputFolder %s'%inputFolder)
#    
#
#
