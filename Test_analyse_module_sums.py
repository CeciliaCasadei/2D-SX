# -*- coding: utf-8 -*-
import pickle
import numpy


results_file_old = open('./Output_imageSums_modules_2/module_results_2_11_totTerms_1353.pkl', 'rb')

orbit_results_old = pickle.load(results_file_old)

print orbit_results_old.shape

nUsed = 0
sigma_Xs = []
sigma_Ys = []
x0s = []
y0s = []
for module in orbit_results_old:
    nTerms = module[6]
    nUsed = nUsed + nTerms
    sigma_Xs.append(abs(float(module[2])))
    sigma_Ys.append(abs(float(module[3])))
    x0s.append(float(module[4]))
    y0s.append(float(module[5]))
    
    
print nUsed
avg_sigma_X = numpy.average(sigma_Xs)
print avg_sigma_X
avg_sigma_Y = numpy.average(sigma_Ys)
print avg_sigma_Y
delta_x = max(x0s) - min(x0s)
print max(x0s), min(x0s), delta_x
delta_y = max(y0s) - min(y0s)
print max(y0s), min(y0s), delta_y

results_file_new = open('./Output_imageSums_modules/module_results_2_11_nUsed_1087_nTot_1342.pkl', 'rb')

orbit_results_new = pickle.load(results_file_new)

print orbit_results_new.shape

nUsed = 0
sigma_Xs = []
sigma_Ys = []
x0s = []
y0s = []
for module in orbit_results_new:
    nTerms = module[6]
    nUsed = nUsed + nTerms
    sigma_Xs.append(abs(float(module[2])))
    sigma_Ys.append(abs(float(module[3])))
    x0s.append(float(module[4]))
    y0s.append(float(module[5]))
    
    
print nUsed
avg_sigma_X = numpy.average(sigma_Xs)
print avg_sigma_X
avg_sigma_Y = numpy.average(sigma_Ys)
print avg_sigma_Y
delta_x = max(x0s) - min(x0s)
print max(x0s), min(x0s), delta_x
delta_y = max(y0s) - min(y0s)
print max(y0s), min(y0s), delta_y
