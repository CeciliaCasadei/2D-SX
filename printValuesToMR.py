# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import joblib

folder = './Output_runMergingVsModel/Shannon_sampling'
FW_data = joblib.load('%s/French_Wilson/FW_uniques.jbl'%folder)
print FW_data.shape


fOpen = open('%s/h_k_l_F_sigF_FrenchWilson.txt'%folder, 'w')


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
