# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import numpy
#from numpy.polynomial import polynomial as P


import braggRodClass
import makeOrbits
#from model_poly_order import poly_order


def mergingFunction(myArguments):
    input_str = '--inputFolder <inputFolder> --resolutionLimit <resolutionLimit>'
    
    # DEFAULTS
    #runNumbers = ['0195', '0196', '0197', '0200', '0201']
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']

    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["inputFolder=", 
                                               "resolutionLimit="])
    except getopt.GetoptError:
        print 'Usage: python merging.py %s'%input_str
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python merging.py %s'%input_str
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
    
    # PREPARE FOLDERS
    outputFolder = '%s/mergedRods'%inputFolder      
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    rodsFolder = '%s/braggRodObjects'%inputFolder
    if not os.path.exists(rodsFolder):
        os.mkdir(rodsFolder)
#    if not os.path.exists('%s/intensitiesAtZero'%inputFolder):
#        os.mkdir('%s/intensitiesAtZero'%inputFolder)
#    if not os.path.exists('%s/model'%inputFolder):
#        os.mkdir('%s/model'%inputFolder)
        
    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
        
    print '%d Rods'%len(rodIndices)   
                  
#    IsAtZero = numpy.zeros(shape=(len(rodIndices), 3))
#    lattice_model = []
    nPos = 0
    nNeg = 0
    rodNumber = 0   
    # PLOT I vs QROD, FOR EVERY ROD, AND PRODUCE MODEL (POLYNOMIAL FIT)    
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector = []
        Irod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I) POINTS FROM ALL RUNS 
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            listDirectory = '%s/spotsMatricesList-Scaled-r%s'%(inputFolder, 
                                                               runNumber)
            myList = joblib.load('%s/r%s_scaledSpotsMatricesList.jbl'
                                  %(listDirectory, runNumber))
        
            for latticeMatrix in myList: # h k qRod I flag
                latticeMatrix = numpy.asarray(latticeMatrix)
                if latticeMatrix[0, 4] == 1:
                    for spot in latticeMatrix:
                        h = spot[0]
                        k = spot[1]
                        if ((h == hRod       and k == kRod) or 
                            (h == -hRod-kRod and k == hRod) or 
                            (h == kRod       and k == -hRod-kRod)):
                            Irod_vector.append(spot[3])
                            Qrod_vector.append(spot[2])
                            if spot[3] >= 0:
                                nPos = nPos+1
                            else:
                                nNeg = nNeg+1
                        if ((h == -hRod     and k == -kRod) or 
                            (h == hRod+kRod and k == -hRod) or 
                            (h == -kRod     and k == hRod+kRod)):
                            Irod_vector.append(spot[3])
                            Qrod_vector.append(-spot[2])
                            if spot[3] >= 0:
                                nPos = nPos+1
                            else:
                                nNeg = nNeg+1
         
        # REMOVE NAN VALUES
        cleanedList_Irod = [Irod_vector[i] for i in range(0, len(Irod_vector)) 
                                           if not numpy.isnan(Irod_vector[i])]
        cleanedList_Qrod = [Qrod_vector[i] for i in range(0, len(Irod_vector)) 
                                           if not numpy.isnan(Irod_vector[i])]
                
#        # POLYNOMIAL FIT ORDER
#        n = int(poly_order['[%s, %s]'%(hRod, kRod)])
#        
#        # POLYNOMIAL FIT
#        c, stats = P.polyfit(cleanedList_Qrod, cleanedList_Irod, n, full=True) 
        Qrod_min = min(cleanedList_Qrod)
        Qrod_max = max(cleanedList_Qrod)
#        x_fit = numpy.linspace(Qrod_min, Qrod_max, 
#                               num=(Qrod_max-Qrod_min)/0.001, 
#                               endpoint = True)
#        y_fit = numpy.zeros(shape=x_fit.shape)
#        for exponent in range(0, n+1):
#            y_fit = y_fit + c[exponent]*(x_fit**exponent)
#            
#        for i in range(0, len(x_fit)):
#            spot_model = [int(hRod), int(kRod), x_fit[i], y_fit[i]] # h k qRod I
#            lattice_model.append(spot_model)
               
        # PLOT BRAGG ROD
        matplotlib.pyplot.scatter(cleanedList_Qrod, 
                                  cleanedList_Irod, 
                                  marker='o', 
                                  color='c', 
                                  alpha = 0.15, 
                                  s=10)
        #matplotlib.pyplot.plot(x_fit, y_fit, '.b-')
        myAxis = matplotlib.pyplot.gca()
        matplotlib.pyplot.axhline(y=0,  linewidth=0.5, color = 'b')
        matplotlib.pyplot.axhline(y=10, linewidth=0.5, color = 'b')
        myAxis.set_xlim([-0.60,+0.60])
        posScale = 1.1*max(cleanedList_Irod)
        negScale = 1.1*min(cleanedList_Irod)
        myAxis.set_ylim([negScale,posScale])
        myAxis.set_xlabel("q$_z$ (A$^{-1}$)", fontsize = 12, 
                                              rotation = 'horizontal')
        matplotlib.pyplot.savefig('%s/mergedRod_%d_%d.png'
                                   %(outputFolder, 
                                     hRod, 
                                     kRod))
        matplotlib.pyplot.close()
        
        # COLLECT I(qRod = 0)
#        IsAtZero[rodNumber, 0]= hRod
#        IsAtZero[rodNumber, 1]= kRod
#        IsAtZero[rodNumber, 2]= c[0]
        
        # GENERATE braggRod OBJECT
        braggRodObject = braggRodClass.braggRod(hRod, 
                                                kRod, 
                                                Qrod_min, 
                                                Qrod_max)
        braggRodObject.setExperimentalPoints(cleanedList_Qrod, 
                                             cleanedList_Irod)
#        braggRodObject.setModelPoints(x_fit, 
#                                      y_fit)
#        braggRodObject.setModelCoefficients(c)
        joblib.dump(braggRodObject, '%s/braggRodObject_%d_%d.jbl'
                                     %(rodsFolder, hRod, kRod))
                                     
        rodNumber = rodNumber + 1
        
#    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
#            
#    joblib.dump(IsAtZero, 
#                '%s/intensitiesAtZero/intensitiesAtZero.jbl'%inputFolder)
#    joblib.dump(lattice_model, 
#                '%s/model/lattice_model.jbl'%inputFolder)
    print 'N_NEG / N_TOT = %d/%d = %.2f'%(nNeg, 
                                          nPos+nNeg, 
                                          float(nNeg)/(nPos+nNeg))

if __name__ == "__main__":
    print "\n**** CALLING merging ****"
    mergingFunction(sys.argv[1:])    