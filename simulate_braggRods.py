# -*- coding: utf-8 -*-
import os
import numpy
import matplotlib.pyplot

import makeOrbits

twoD_resolution = 3.5
folderName = './braggRodsSimulations'
#outputFolder = './%s/Figures_SFALL'%folderName
#outputFolder = './%s/Figures_PHENIX'%folderName
outputFolder = './%s/Figures_SFALL_test_movingChain'%folderName
c = 504.5 #A

c_star = 2*numpy.pi/c


if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)
orbits = makeOrbits.makeOrbitsFunction(twoD_resolution)

N = 0
sum_of_ratios = 0
total_I_1fbb = 0
total_I_1fbk = 0

for orbit in orbits:
    h_label = orbit.label[0]
    k_label = orbit.label[1]
    
    if h_label >= 0 and k_label >= 0:
        print h_label, k_label
        
        matplotlib.pyplot.figure()        
        
        #SF_1fbb_open = open('%s/SF_from_1fbb_c_504p5A_ed.txt'%folderName, 'r')  #SFALL PDB -> Fs
        #SF_1fbb_open = open('%s/1fbb_edited_fmodel_ed.txt'%folderName, 'r')     #PHENIX PDB -> Fs
        SF_1fbb_open = open('%s/SF_from_1fbb_c_504p5A_Test.txt'%folderName, 'r') #SFALL TEST MOVING CHAIN
        rod_1fbb_qRod = []
        rod_1fbb_I = []
        for line in SF_1fbb_open:
            h = int(line[0:5])
            k = int(line[5:9])
            if h == h_label and k == k_label:
                l = int(line[9:14])
                qRod = l*c_star
                rod_1fbb_qRod.append(qRod)
                F = float(line[14:25])
                I = F**2/124651
                rod_1fbb_I.append(I)
        matplotlib.pyplot.plot(rod_1fbb_qRod, rod_1fbb_I, 'b')
        SF_1fbb_open.close()
        
        #SF_1fbk_open = open('%s/SF_from_1fbk_c_504p5A_ed.txt'%folderName, 'r')                  #SFALL PDB -> Fs
        #SF_1fbk_open = open('%s/1fbk_edited_fmodel_ed.txt'%folderName, 'r')                     #PHENIX PDB -> Fs
        SF_1fbk_open = open('%s/SF_from_1fbb-modified_chains_c_504p5A_Test.txt'%folderName, 'r') #SFALL TEST MOVING CHAIN
        rod_1fbk_qRod = []
        rod_1fbk_I = []
        for line in SF_1fbk_open:
            h = int(line[0:5])
            k = int(line[5:9])
            if h == h_label and k == k_label:
                l = int(line[9:14])
                qRod = l*c_star
                rod_1fbk_qRod.append(qRod)
                F = float(line[14:25])
                I = F**2/124486
                rod_1fbk_I.append(I)
        matplotlib.pyplot.plot(rod_1fbk_qRod, rod_1fbk_I, 'r')
        SF_1fbk_open.close()
        
        for i in range(0, len(rod_1fbb_qRod)):
            qRod_1fbb = rod_1fbb_qRod[i]
            qRod_1fbk = rod_1fbk_qRod[i]
            if qRod_1fbb != qRod_1fbk:
                print 'Diff'
                
            I_1fbb = rod_1fbb_I[i]
            total_I_1fbb = total_I_1fbb + I_1fbb
            I_1fbk = rod_1fbk_I[i]
            total_I_1fbk = total_I_1fbk + I_1fbk
            delta_I = abs(I_1fbb - I_1fbk)
            avg_I = (I_1fbb+I_1fbk)/2
            ratio = delta_I/abs(avg_I)
            sum_of_ratios = sum_of_ratios + ratio
            N = N+1
            
        matplotlib.pyplot.savefig('%s/rod_%d_%d.png'%(outputFolder, h_label, k_label), dpi=2*96)
        matplotlib.pyplot.close()
avg_delta = sum_of_ratios/N
        
print avg_delta, N
print total_I_1fbb/N, total_I_1fbk/N