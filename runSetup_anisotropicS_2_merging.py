# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_merging.py
    
    
# DETERMINE AND APPLY RUN - RUN TRANSFORMATIONS
deltaQrodThreshold = 0.005   # A-1
n_minThreshold = 6
nLatticePairs = 1000

flag = 0
if flag == 1:
    os.system('python transform_CCmethod_mergeRuns.py --dQrod %f \
                                                      --nMin %d \
                                                      --nLatticePairs %d'
                                                      %(deltaQrodThreshold, 
                                                        n_minThreshold, 
                                                        nLatticePairs))    
    

### DEPRECATED ###
# DETERMINE AND APPLY RUN - RUN SCALE FACTORS
#deltaQrodThreshold = 0.003    # A-1
#n_minThreshold = 8
#nLatticePairs = 1400
#resolution_3D = 12.0          # A
#
#flag = 0
#if flag == 1:
#    os.system('python scaling_mergeRuns.py --dQrod %f \
#                                           --nMin %d \
#                                           --nLatticePairs %d \
#                                           --resolution_3D %f'
#                                           %(deltaQrodThreshold, 
#                                             n_minThreshold, 
#                                             nLatticePairs, 
#                                             resolution_3D))
                                             
                                             

resolution_3D = 18.0  # A
cellSize = 62.45      # A

flag = 0
if flag == 1:
    os.system('python scaling_mergeRuns_total_I.py --resolution_3D %f \
                                                   --cellSize %f'
                                                   %(resolution_3D,
                                                     cellSize))
    
    
    
# PLOT MERGED RODS, BUILD SINC MODEL
resolutionLimit = 6.0      # A. 2D
thickness = 45             # A
damping = 80
folder = './Output_runMerging'

flag = 0
if flag == 1:
    os.system('python merging.py --inputFolder %s \
                                 --resolutionLimit %f'
                                 %(folder, 
                                   resolutionLimit))
                                   
                                   
                                   
flag = 1
if flag == 1:
    os.system('python rodsFit_shannonTheo.py \
               --resolutionLimit %f \
               --thickness %f \
               --damping %f \
               --folder %s'
               %(resolutionLimit, 
                 thickness, 
                 damping, 
                 folder))



newFolder = './Output_runMergingVsModel'
if not os.path.exists(newFolder):
    os.mkdir(newFolder)   