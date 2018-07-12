# -*- coding: utf-8 -*-
import joblib

folder = './Output_runMergingVsModel/Shannon_sampling'
R_data = joblib.load('%s/Rfactor/R_factor_bins.jbl'%folder)
CC_data = joblib.load('%s/CChalf/CChalf_bins.jbl'%folder)
StoN_data = joblib.load('%s/French_Wilson/FW_bins.jbl'%folder)
print R_data.shape, CC_data.shape, StoN_data.shape


fOpen = open('%s/stats.txt'%folder, 'w')
fOpen.write('Res 3D: low-high    Nobs    Nunique Redundancy      R  CChalf  CCstar (S/N)_FW\n')
nBins = R_data.shape[0]
for bin_line in range(0, nBins-1):
    
    low  = R_data[bin_line][0]
    high = R_data[bin_line][1]
    Nobs = R_data[bin_line][2]
    R    = R_data[bin_line][3]
    
    Nunique = CC_data[bin_line][2]
    CChalf  = CC_data[bin_line][3]
    CCstar  = CC_data[bin_line][4]
    
    StoN = StoN_data[bin_line][3]
    fOpen.write('%6.2f & %6.2f & %10d & %10d & %10.2f & %6.2f & %6.4f & %6.4f & %6.2f \\ \n'%(low,
                                                                         high,
                                                                         Nobs,
                                                                         Nunique,
                                                                         (float(Nobs))/Nunique,
                                                                         R,
                                                                         CChalf,
                                                                         CCstar,
                                                                         StoN))
    
    
fOpen.close()