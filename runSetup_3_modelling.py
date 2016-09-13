# -*- coding: utf-8 -*-
import os
import sys


# SETUP FOR CURRENT RUN                                                                                                                                                                       
if len(sys.argv) != 2:
    print("[USAGE] %s runnumber" % sys.argv[0])
    sys.exit(-1)

runNumber = '%04d' % int(sys.argv[1])



# DETERMINE LATTICES TO MODEL TRANSFORMATIONS
deltaQrodThreshold = 0.001
n_minThreshold = 6
nUsedLattices = 'all'
nTriangles = 100
nGoodFraction = 0.7

flag = 1
if flag == 1:
    os.system('python model_transformVsModel.py --runNumber %s --dQrod %f --nMin %d --nLattices %s --nTriangles %d --nGoodFraction %f'
    %(runNumber, deltaQrodThreshold, n_minThreshold, nUsedLattices, nTriangles, nGoodFraction))
    

    
# APPLY INDICES TRANSFORMATIONS
flag = 1
if flag == 1:
    os.system('python model_applyTransformations.py --runNumber %s'%runNumber)



# DETERMINE LATTICE TO MODEL SCALE FACTOR    
flag = 1
if flag == 1:
    os.system('python model_scaleVsModel.py --runNumber %s'%runNumber)