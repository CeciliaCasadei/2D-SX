# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot
#PDB_code = '2ntu'
PDB_code = 'SF_from_MR_model'
SF_file = open('../PDB_models/%s_only_SF.txt'%PDB_code, 'r') # File generated by SFALL: h k l f phi(F)
                                                             # Starting from PDB file, after editing UC, SG and SCALE matrix
c = 100.9 # A
c_long = 5*c
c_long_star = 2*numpy.pi/c_long

extractedReflections = []
for reflection in SF_file:
    split_reflection = reflection.split()
    h = int(split_reflection[0])
    k = int(split_reflection[1])
    l = int(split_reflection[2])
    F = float(split_reflection[3])
    I = F**2
    qRod = l*c_long_star
    print '%d %d %d %.2f %.2f %.2f'%(h, k, l, F, qRod, I)
    extractedReflection = [h, k, l, F, qRod, I]
    extractedReflections.append(extractedReflection)
    
extractedReflections = numpy.asarray(extractedReflections)

rodIndices = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
                  
for rodIndices_singleRod in rodIndices:
    F_rod = []
    I_rod = []
    q_rod = []
    hRod = rodIndices_singleRod[0]
    kRod = rodIndices_singleRod[1]
    print '%d %d'%(hRod, kRod)
    for extractedReflection in extractedReflections:
        if extractedReflection[0] == hRod and extractedReflection[1] == kRod:
            I_rod.append(extractedReflection[5])
            q_rod.append(extractedReflection[4])
            F_rod.append(extractedReflection[3])
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(q_rod, I_rod, '-o')
    matplotlib.pyplot.savefig('../PDB_models/%s_I_rod_%d_%d.png'%(PDB_code, hRod, kRod))
    matplotlib.pyplot.close()    
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(q_rod, F_rod, '-o')
    myAxis = matplotlib.pyplot.gca()
    myAxis.set_xlim([-0.32, +0.32])
    matplotlib.pyplot.savefig('../PDB_models/%s_F_rod_%d_%d.png'%(PDB_code, hRod, kRod))
    matplotlib.pyplot.close() 
    