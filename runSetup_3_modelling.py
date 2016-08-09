# -*- coding: utf-8 -*-
import os
import sys


# SETUP FOR CURRENT RUN                                                                                                                                                                       
if len(sys.argv) != 2:
    print("[USAGE] %s runnumber" % sys.argv[0])
    sys.exit(-1)

runNumber = '%04d' % int(sys.argv[1])

# DETERMINE LATTICES TO MODEL TRANSFORMATIONS

flag = 1
if flag == 1:
    os.system('python model_transformVsModel.py --runNumber %s'%runNumber)
