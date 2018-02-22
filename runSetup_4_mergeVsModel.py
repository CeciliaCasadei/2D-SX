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
    os.system('python merging.py --inputFolder %s --resolutionLimit %f'
              %(inputFolder, resolutionLimit))
              
              
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
  

             
# CALCULATE RESOLUTION LIMITS AND BINS          
resolutionLimit_2D = 6.0
inputFolder = './Output_runMergingVsModel/Shannon_sampling'
cellSize = 62.45

flag = 0
if flag == 1:
    os.system('python dataLimits.py --inputFolder %s \
                                    --resolutionLimit %f \
                                    --cellSize %f '
                                    %(inputFolder, 
                                      resolutionLimit_2D,
                                      cellSize))
               
# CALCULATE R-FACTOR
flag = 0
if flag == 1:
    os.system('python calculate_Rfactor.py --inputFolder %s \
                                           --thickness %f \
                                           --damping %f \
                                           --cellSize %f'
                                           %(inputFolder, 
                                             thickness,
                                             damping,
                                             cellSize))
                                             
                                             
                                             
# CALCULATE CChalf
flag = 0
if flag == 1:
    os.system('python calculate_CChalf.py --inputFolder %s \
                                          --thickness %f \
                                          --damping %f \
                                          --cellSize %f'
                                          %(inputFolder, 
                                            thickness,
                                            damping,
                                            cellSize))
                                            
 

# FRENCH-WILSON METHOD
flag = 0
if flag == 1:
    os.system('python calculate_FW_values.py --inputFolder %s \
                                             --thickness %f \
                                             --damping %f \
                                             --cellSize %f'
                                             %(inputFolder, 
                                               thickness,
                                               damping,
                                               cellSize)) 
                                               
                                               
                                               
# PRINT ALL STATS
flag = 0
if flag == 1:
    os.system('python printStats.py')
    
    
                                               
# FROM FW RESULTS, EXTRACT DATA TO MR  
flag = 0
if flag == 1:
    os.system('python printValuesToMR.py')
    
    
    
# EXTRACT DATA TO MR, FROM SHANNON FITS, NO FRENCH-WILSON
flag = 0
if flag == 1:
    os.system('python prepare_MRdata.py --resolutionLimit %f \
                                        --thickness %f \
                                        --damping %f'
                                        %(resolutionLimit,
                                          thickness,
                                          damping))
            
                                                                                                                   
                                            
# MODEL COMPLETENESS VS TILT ANGLE
flag = 0
if flag == 1:
    os.system('python simulate_completeness.py')