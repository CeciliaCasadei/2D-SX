# -*- coding: utf-8 -*-
#import numpy
#import matplotlib.pyplot
#
#import makeOrbits
#
#def Correlate(x1, x2):
#    x1Avg = numpy.average(x1)
#    x2Avg = numpy.average(x2)
#    numTerm = numpy.multiply(x1-x1Avg, x2-x2Avg)
#    num = numTerm.sum()
#    resX1Sq = numpy.multiply(x1-x1Avg, x1-x1Avg)
#    resX2Sq = numpy.multiply(x2-x2Avg, x2-x2Avg)
#    den = numpy.sqrt(numpy.multiply(resX1Sq.sum(), resX2Sq.sum()))
#    CC = num/den
#    return CC

### CHECK LATTICE TRANSFORMATIONS ###
runNumber = '0127'
rigid_transformation = 'ip'
Ts_py = []
Ts_EMC = []
transformationsFile_python = open('./Cmp_EMC_py/r%s-finalOrientations.txt'%runNumber, 'r')    
transformationsFile_python_list = list(transformationsFile_python)
transformationsFile_EMC = open('./Cmp_EMC_py/model10_r%s_sym.txt'%runNumber, 'r')                                    
transformationsFile_EMC_list = list(transformationsFile_EMC)
n_good = 0
n_bad = 0
n_nan = 0
for lineNumber in range(0, len(transformationsFile_EMC_list)):
    T_line_EMC = transformationsFile_EMC_list[lineNumber]    
    T_line_EMC_split = T_line_EMC.split(',')
    T_EMC = int(T_line_EMC_split[1]) 
    
    T_line_py = transformationsFile_python_list[lineNumber]    
    T_line_py_split = T_line_py.split()
    T_py = T_line_py_split[4]
    
    
    ### ip between py and EMC ##
    if rigid_transformation == 'ip':
        try:
            T_py = int(T_py)
            if T_py == 1:                                   ### i
                if T_EMC == 6 or T_EMC == 8 or T_EMC == 10: ### p
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 0:                                   ### I
                if T_EMC == 7 or T_EMC == 9 or T_EMC == 11: ### ip
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 2:                                   ### p
                if T_EMC == 1 or T_EMC == 3 or T_EMC == 5:  ### i
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 3:                                   ### ip
                if T_EMC == 0 or T_EMC == 4 or T_EMC == 2:  ### I
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
        except:
            n_nan = n_nan + 1
            
    ### I between py and EMC ##
    if rigid_transformation == 'I':
        try:
            T_py = int(T_py)
            if T_py == 2:                                   ### p
                if T_EMC == 6 or T_EMC == 8 or T_EMC == 10: ### p
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 3:                                   ### ip
                if T_EMC == 7 or T_EMC == 9 or T_EMC == 11: ### ip
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 1:                                   ### i
                if T_EMC == 1 or T_EMC == 3 or T_EMC == 5:  ### i
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 0:                                   ### I
                if T_EMC == 0 or T_EMC == 4 or T_EMC == 2:  ### I
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
        except:
            n_nan = n_nan + 1
            
    ### i between py and EMC ##
    if rigid_transformation == 'i':
        try:
            T_py = int(T_py)
            if T_py == 3:                                   ### ip
                if T_EMC == 6 or T_EMC == 8 or T_EMC == 10: ### p
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 2:                                   ### p
                if T_EMC == 7 or T_EMC == 9 or T_EMC == 11: ### ip
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 0:                                   ### I
                if T_EMC == 1 or T_EMC == 3 or T_EMC == 5:  ### i
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 1:                                   ### i
                if T_EMC == 0 or T_EMC == 4 or T_EMC == 2:  ### I
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
        except:
            n_nan = n_nan + 1
            
            
    ### p between py and EMC ##
    if rigid_transformation == 'p':
        try:
            T_py = int(T_py)
            if T_py == 0:                                   ### I
                if T_EMC == 6 or T_EMC == 8 or T_EMC == 10: ### p
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 1:                                   ### i
                if T_EMC == 7 or T_EMC == 9 or T_EMC == 11: ### ip
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 3:                                   ### ip
                if T_EMC == 1 or T_EMC == 3 or T_EMC == 5:  ### i
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
            if T_py == 2:                                   ### p
                if T_EMC == 0 or T_EMC == 4 or T_EMC == 2:  ### I
                    n_good = n_good + 1
                else:
                    n_bad = n_bad + 1
        except:
            n_nan = n_nan + 1
        
print 'GOOD: %d' %n_good
print 'BAD: %d' %n_bad
print 'NAN: %d' %n_nan


### RUN 0127 ###
#resolutionLimit = 7.0    
#
#Is_python = []
#Is_EMC = []
#
#intensityFile_python = open('./Cmp_EMC_py/mergedIntensities_h_k_qRod_l_avgI_sigmaI.txt', 'r')       # h k qRod l I sig(I)  - P3 merged   - to 4A
#intensityFile_python_list = list(intensityFile_python)
#intensityFile_EMC = open('./Cmp_EMC_py/model10_non_zero_l_sym.hkl', 'r')                                        # h k l I sig(I) nmeas - P3 unmerged - to 7A
#intensityFile_EMC_list = list(intensityFile_EMC)
#
## DEFINE ROD INDICES       
#orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
#
#for orbit in orbits:               # Each orbit contains 3 P3 related reflections
#    orbit_label = orbit.label
#    h_orbit = orbit_label[0]
#    k_orbit = orbit_label[1]
#    print h_orbit, k_orbit   
#    
#    # EXTRACT PYTHON INTENSITY
#    for I_line_py in intensityFile_python_list:
#        splittedLine = I_line_py.split()
#        h = int(splittedLine[0])
#        k = int(splittedLine[1])
#        I_python = float(splittedLine[4])
#        if k_orbit != 0:
#            if h == -k_orbit and k == -h_orbit:  ### TRANSFORMATION ip PYTHON - EMC ###
#                Is_python.append(I_python)
#                print h_orbit, k_orbit, I_python
#        else:
#            if h == h_orbit and k == k_orbit:    ### TRANSFORMATION ip PYTHON - EMC ###
#                Is_python.append(I_python)
#                print h_orbit, k_orbit, I_python
#    
#    # EXTRACT EMC INTENSITY           
#    I_EMC_unmerged = []            
#    for I_line_EMC in intensityFile_EMC_list:
#        splittedLine = I_line_EMC.split()
#        h = int(splittedLine[0])
#        k = int(splittedLine[1])
#        l = int(splittedLine[2])
#        I_EMC = float(splittedLine[3])
#        if l < 0:
#            h = -h
#            k = -k
#            l = -l
#        if (h == h_orbit and k == k_orbit) or (h == -h_orbit-k_orbit and k == h_orbit) or (h == k_orbit and k == -h_orbit-k_orbit):
#            I_EMC_unmerged.append(I_EMC)
#    I_EMC_unmerged = numpy.asarray(I_EMC_unmerged)
#    print I_EMC_unmerged
#    I_EMC_merged = numpy.average(I_EMC_unmerged)
#    Is_EMC.append(I_EMC_merged)
#    print I_EMC_merged
#        
#        
#print len(Is_python), len(Is_EMC)
#
#CC = Correlate(Is_python, Is_EMC)
#print 'CC python-EMC: %.5f'%CC
#
#intensityFile_EMC.close()
#intensityFile_python.close()
#
#
#matplotlib.pyplot.figure()
#matplotlib.pyplot.title('n pairs = %d - CC = %.8f'%(len(Is_python), CC), y=1.05)
#matplotlib.pyplot.scatter(Is_python, Is_EMC)
#matplotlib.pyplot.gca().set_xscale('log')
#matplotlib.pyplot.gca().set_yscale('log')
#matplotlib.pyplot.gca().set_xlim([0.1, 1000])
#matplotlib.pyplot.gca().set_ylim([0.1, 1000])
#matplotlib.pyplot.gca().set_ylabel(r'I - EMC')
#matplotlib.pyplot.gca().set_xlabel(r'I - PYTHON')
#matplotlib.pyplot.axes().set_aspect('equal')
#matplotlib.pyplot.savefig('./Cmp_EMC_py/r0127_Is_python_VS_Is_EMC_7A_logscale_non_zero_l_sym.png', dpi=96*4)
#matplotlib.pyplot.close()
