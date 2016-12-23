# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot



def Correlate(x1, x2):
    x1Avg = numpy.average(x1)
    x2Avg = numpy.average(x2)
    numTerm = numpy.multiply(x1-x1Avg, x2-x2Avg)
    num = numTerm.sum()
    resX1Sq = numpy.multiply(x1-x1Avg, x1-x1Avg)
    resX2Sq = numpy.multiply(x2-x2Avg, x2-x2Avg)
    den = numpy.sqrt(numpy.multiply(resX1Sq.sum(), resX2Sq.sum()))
    CC = num/den
    return CC
    




Fs_obs_file_python = open('./Cmp_Fo_Fc/experimentalReflections_c_504p5A_F_sigF.hkl', 'r')       # h k l F sig(F)  - P3 merged   - to 4A
Fs_obs_file_python_list = list(Fs_obs_file_python)
Fs_calc_file_phenix = open('./Cmp_Fo_Fc/MR_bR.1_edit_refine_001_refine_003_f_model_ed.txt', 'r')                                        
Fs_calc_file_phenix_list = list(Fs_calc_file_phenix)
Fs_calFs_calc_file_phenix_write = open('./Cmp_Fo_Fc/MR_bR.1_edit_refine_001_refine_003_f_model_ed_short.txt', 'w')       

print len(Fs_obs_file_python_list)
print len(Fs_calc_file_phenix_list)
print len(Fs_calc_file_phenix_list)/5

for Fcalc_line in Fs_calc_file_phenix_list:
    if len(Fcalc_line.split()) == 9 :
        Fs_calFs_calc_file_phenix_write.write(Fcalc_line)  # h k l Fo sig(Fo) RfreeFlag Fmodel PHImodel Fcalc
        
        
Fs_obs_file_python.close()
Fs_calc_file_phenix.close()
Fs_calFs_calc_file_phenix_write.close()

Fs_calc_file_phenix_short = open('./Cmp_Fo_Fc/MR_bR.1_edit_refine_001_refine_003_f_model_ed_short.txt', 'r')       
Fs_calc_file_phenix_short_list = list(Fs_calc_file_phenix_short)
Fs_calc_file_phenix_short.close()
print len(Fs_calc_file_phenix_short_list)
        

# CMP phenix Fobs - phenix Fmodel
Fs_obs_phenix = []
Fs_model_phenix = []
Fs_calc_phenix = []
for Fcalc_line in Fs_calc_file_phenix_short_list:
    F_obs_phenix = float(Fcalc_line.split()[3])
    Fs_obs_phenix.append(F_obs_phenix)
    F_model_phenix = float(Fcalc_line.split()[6])
    Fs_model_phenix.append(F_model_phenix)
    F_calc_phenix = float(Fcalc_line.split()[8])
    Fs_calc_phenix.append(F_calc_phenix)
    
print len(Fs_obs_phenix), len(Fs_model_phenix), len(Fs_calc_phenix)
cc = Correlate(Fs_obs_phenix, Fs_model_phenix)
print 'CC Phenix Fobs - Phenix Fmodel: %f'%cc
cc = Correlate(Fs_obs_phenix, Fs_calc_phenix)
print 'CC Phenix Fobs - Phenix Fcalc: %f'%cc

Fs_obs_python = []
Fs_obs_phenix = []
Fs_model_phenix = []
Fs_calc_phenix = []

# CMP python Fobs - phenix Fobs
for Fobs_line in Fs_obs_file_python_list:
    h = int(Fobs_line.split()[0])
    k = int(Fobs_line.split()[1])
    l = int(Fobs_line.split()[2])
    F_obs_python = float(Fobs_line.split()[3])
    for Fcalc_line in Fs_calc_file_phenix_short_list:
        if int(Fcalc_line.split()[0]) == h and  int(Fcalc_line.split()[1]) == k and int(Fcalc_line.split()[2]) == l:           
            F_obs_phenix = float(Fcalc_line.split()[3])
            F_model_phenix = float(Fcalc_line.split()[6])
            F_calc_phenix = float(Fcalc_line.split()[8])
            
            Fs_calc_phenix.append(F_calc_phenix)
            Fs_model_phenix.append(F_model_phenix)
            Fs_obs_phenix.append(F_obs_phenix)
            Fs_obs_python.append(F_obs_python)
        
print len(Fs_obs_phenix), len(Fs_obs_python), len(Fs_calc_phenix), len(Fs_model_phenix)
cc = Correlate(Fs_obs_phenix, Fs_obs_python)
print 'CC Phenix Fobs - Python Fobs: %f'%cc
cc = Correlate(Fs_calc_phenix, Fs_obs_python)
print 'CC Phenix Fcalc - Python Fobs: %f'%cc 
cc = Correlate(Fs_model_phenix, Fs_obs_python)
print 'CC Phenix Fmodel - Python Fobs: %f'%cc 