# -*- coding: utf-8 -*-
import os
import joblib

folder = './Output_runMergingVsModel/Shannon_sampling'

nBins = 15
if os.path.exists('%s/h_k_l_F_sigF_FrenchWilson.txt'%folder):
    os.remove('%s/h_k_l_F_sigF_FrenchWilson.txt'%folder)
fOpen = open('%s/h_k_l_F_sigF_FrenchWilson.txt'%folder, 'a')

for i in range(0, nBins):
    FW_data = joblib.load('%s/French_Wilson/FW_uniques_bin_%d.jbl'%(folder, i))
    print FW_data.shape   
       
    for spot in FW_data:
        
        h    = int(spot[0])
        k    = int(spot[1])
        l    = int(spot[2])
        F    = spot[8]
        sigF = spot[9]
        
        fOpen.write('%5d %5d %5d %8.4f %8.4f\n'%(h,
                                                 k,
                                                 l,
                                                 F,
                                                 sigF))
    
    
fOpen.close()
