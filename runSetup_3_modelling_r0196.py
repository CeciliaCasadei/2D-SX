# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_3_modelling_rXXXX.py
    
    
# DETERMINE LATTICES TO MODEL TRANSFORMATIONS
runNumber = '0196'

deltaQrodThreshold = 0.001
n_minThreshold = 6
nUsedLattices = 'all'
nTriangles = 100
nGoodFraction = 0.7


flag = 0
if flag == 1:
    os.system('python model_transformVsModel.py --runNumber %s --dQrod %f --nMin %d --nLattices %s --nTriangles %d --nGoodFraction %f'
    %(runNumber, deltaQrodThreshold, n_minThreshold, nUsedLattices, nTriangles, nGoodFraction))
    
    
# APPLY INDICES TRANSFORMATIONS
flag = 0
if flag == 1:
    os.system('python model_applyTransformations.py --runNumber %s'%runNumber)

# DETERMINE LATTICE TO MODEL SCALE FACTOR    
flag = 1
if flag == 1:
    os.system('python model_scaleVsModel.py --runNumber %s'%runNumber)

