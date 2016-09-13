# -*- coding: utf-8 -*-
import joblib
import matplotlib.pyplot
import os
import numpy


if not os.path.exists('./Output_polynomialModel'):
    os.mkdir('./Output_polynomialModel')

hkl_file = open('./Output_polynomialModel/experimentalReflections_c_504p5A_F_sigF.hkl', 'w')

c = 100.9
c_long = 5*c
c_long_star = 2*numpy.pi/c_long


rodIndices_all = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                  
model = joblib.load('./Output_runMerging_bkp/Output_runMergingVsModel/model/lattice_model.jbl')
print model.shape
for rodIndices in rodIndices_all:
    print rodIndices
    hRod = rodIndices[0]
    kRod = rodIndices[1]
    qRod_model = []
    Irod_model = []
    for model_reflection in model:
        if model_reflection[0] == hRod and model_reflection[1] == kRod:
            qRod_model.append(model_reflection[2])
            Irod_model.append(model_reflection[3])
            
    matplotlib.pyplot.plot(qRod_model, Irod_model, '.b-')
    myAxis = matplotlib.pyplot.gca()
    matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axhline(y=10, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
    myAxis.set_xlim([-0.45,+0.45])
    scale = 1.1*max(Irod_model)
    myAxis.set_ylim([-0.1*scale,1*scale])
    myAxis.set_xlabel("q$_z$ (A$^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.savefig('./Output_polynomialModel/polyFit_rod_%d_%d.png'%(hRod, kRod))
    matplotlib.pyplot.close()
    
    qRod_max = max(qRod_model)
    qRod_min = min(qRod_model)
    l_max = int(qRod_max/c_long_star)
    qRod_model = numpy.asarray(qRod_model)
    for l in range(-l_max, l_max+1):
        print l
        qRod = l*c_long_star
        print qRod
        differenceVector = qRod_model-qRod
        differenceVector = abs(differenceVector)
        index = numpy.argmin(differenceVector)
        I = Irod_model[index]
        if I > 0:
            #sigI = numpy.sqrt(I)
            #hkl_file.write('%3d %3d %3d %8.2f %8.2f\n'%(hRod, kRod, l, I, sigI))
            F = numpy.sqrt(I)
            sigI = F
            sigF = numpy.sqrt(I+F)-F
            hkl_file.write('%3d %3d %3d %8.2f %8.2f\n'%(hRod, kRod, l, F, sigF))
            #hkl_file.write('%3d %3d %3d %8.2f\n'%(hRod, kRod, l, sigI))
        print I
        print '\n'
        
hkl_file.close()

# 1469 reflections out of 4210
    
            
            

