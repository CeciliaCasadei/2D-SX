# -*- coding: utf-8 -*-
import os
import sys

from runSetup_settings import nGoodFractions

# SETUP FOR CURRENT RUN                                                                                                                                                                       
if len(sys.argv) != 2:
    print("[USAGE] %s runnumber" % sys.argv[0])
    sys.exit(-1)

runNumber = '%04d' % int(sys.argv[1])



# DETERMINE LATTICES TO MODEL TRANSFORMATIONS
model = './Output_runMerging/Shannon_sampling/model/lattice_model.jbl'
deltaQrodThreshold = 0.001
n_minThreshold = 6
nUsedLattices = 'all'
nTriangles = 100
nGoodFraction = float(nGoodFractions['%s'%runNumber])  

flag = 0
if flag == 1:
    os.system('python model_transformVsModel.py --runNumber %s \
                                                --model %s \
                                                --dQrod %f \
                                                --nMin %d \
                                                --nLattices %s \
                                                --nTriangles %d \
                                                --nGoodFraction %f'
                                                %(runNumber, 
                                                  model, 
                                                  deltaQrodThreshold, 
                                                  n_minThreshold, 
                                                  nUsedLattices, 
                                                  nTriangles, 
                                                  nGoodFraction))
    

    
# APPLY INDICES TRANSFORMATIONS
flag = 0
if flag == 1:
    os.system('python model_applyTransformations.py --runNumber %s'%runNumber)



# DETERMINE LATTICE TO MODEL SCALE FACTOR 
resolution_3D = 7.0 # A
CC_threshold = 0.92
   
flag = 1
if flag == 1:
    os.system('python model_scaleVsModel.py --runNumber %s \
                                            --resolution_3D %f \
                                            --CC_threshold %f \
                                            --model %s'
                                            %(runNumber, 
                                              resolution_3D,
                                              CC_threshold,
                                              model))