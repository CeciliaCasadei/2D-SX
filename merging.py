# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import numpy
from numpy.polynomial import polynomial as P

import braggRodClass
import makeOrbits

def mergingFunction(myArguments):
    
    # DEFAULTS
    inputFolder = './Output_runMerging'
    runNumbers = ['0127', '0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    
#    rodIndices = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
#                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
#                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder=", "resolutionLimit="])
    except getopt.GetoptError:
        print 'Usage: python merging.py --inputFolder <inputFolder> --resolutionLimit <resolutionLimit>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python merging.py --inputFolder <inputFolder> --resolutionLimit <resolutionLimit>'
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
    
    # PREPARE FOLDERS
    outputFolder = '%s/mergedRods'%inputFolder      
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
    if not os.path.exists('%s/braggRodObjects'%inputFolder):
        os.mkdir('%s/braggRodObjects'%inputFolder)
    if not os.path.exists('%s/intensitiesAtZero'%inputFolder):
        os.mkdir('%s/intensitiesAtZero'%inputFolder)
    if not os.path.exists('%s/model'%inputFolder):
        os.mkdir('%s/model'%inputFolder)
        
    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
        
    print '%d Rods'%len(rodIndices)   
                  
    IsAtZero = numpy.zeros(shape=(len(rodIndices), 3))
    lattice_model = []
    rodNumber = 0   
    # PLOT I vs QROD, FOR EVERY ROD, AND PRODUCE MODEL (POLYNOMIAL FIT OF MEDIAN VALUES)    
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector = []
        Irod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I) POINTS FROM ALL RUNS (AFTER RUN SCALING, WITH AVG SCALE SET TO 1)              
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(inputFolder, runNumber, runNumber))
        
            for latticeMatrix in myList:   # h k qRod I flag
                latticeMatrix = numpy.asarray(latticeMatrix)
                if latticeMatrix[0, 4] == 1:
                    for spot in latticeMatrix:
                        h = spot[0]
                        k = spot[1]
                        if (h == hRod and k == kRod) or (h == -hRod-kRod and k == hRod) or (h == kRod and k == -hRod-kRod):
                            Irod_vector.append(spot[3])
                            Qrod_vector.append(spot[2])
                        if (h == -hRod and k == -kRod) or (h == hRod+kRod and k == -hRod) or (h == -kRod and k == hRod+kRod):
                            Irod_vector.append(spot[3])
                            Qrod_vector.append(-spot[2])
         
        # REMOVE NAN VALUES
        cleanedList_Irod = [Irod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
        cleanedList_Qrod = [Qrod_vector[i] for i in range(0, len(Irod_vector)) if not numpy.isnan(Irod_vector[i])]
                    
        # BINNING
        Qrod_min = min(cleanedList_Qrod)
        Qrod_max = max(cleanedList_Qrod)
        bins, step = numpy.linspace(Qrod_min, Qrod_max, num=(Qrod_max-Qrod_min)/0.008, endpoint = True, retstep = True)
        
        Xs = []
        Ys_means = []
        Ys_medians = []
        # FOR EACH BIN, COLLECT BIN CENTER ON HORIZONTAL AXIS, MEAN AND MEDIAN OF INTENSITY. CHECK WHETHER THE I DISTRIBUTION IN THE BIN IS POISSONIAN.
        for i in range(0, len(bins)-1):
            edge_l = bins[i]
            edge_r = bins[i+1]
            X = (edge_r + edge_l)/2  # Bin center

            binList_Irod = [cleanedList_Irod[item] for item in range(0, len(cleanedList_Irod)) if edge_l <= cleanedList_Qrod[item] <= edge_r]
            
            if len(binList_Irod) > 0:
                Xs.append(X)
                Y_mean = numpy.average(binList_Irod)
                Ys_means.append(Y_mean)
                Y_median = numpy.median(binList_Irod)
                Ys_medians.append(Y_median)
                print '\nN of points in the bin = %d, Mean = %.2f, Median = %.2f'%(len(binList_Irod), Y_mean, Y_median)
                
        
        # POLYNOMIAL FIT ORDER
        if indices == [1, 0]:
            n = 6  
        if indices == [1, 1]:
            n = 6                 
        if indices == [1, 2]:
            n = 5                
        if indices == [1, 3]:
            n = 7              
        if indices == [1, 4]:
            n = 10               
        if indices == [1, 5]:
            n = 13
        if indices == [1, 6]:
            n = 14
        if indices == [1, 7]:
            n = 13
        if indices == [2, 0]:
            n = 4
        if indices == [2, 1]:
            n = 5
        if indices == [2, 2]:
            n = 6
        if indices == [2, 3]:
            n = 13
        if indices == [2, 4]:
            n = 12
        if indices == [2, 5]:
            n = 16
        if indices == [2, 6]:
            n = 12
        if indices == [3, 0]:
            n = 6
        if indices == [3, 1]:
            n = 7
        if indices == [3, 2]:
            n = 13
        if indices == [3, 3]:
            n = 12
        if indices == [3, 4]:
            n = 9
        if indices == [3, 5]:
            n = 10
        if indices == [4, 0]:
            n = 6
        if indices == [4, 1]:
            n = 10
        if indices == [4, 2]:
            n = 12
        if indices == [4, 3]:
            n = 9
        if indices == [4, 4]:
            n = 15
        if indices == [5, 0]:
            n = 14
        if indices == [5, 1]:
            n = 13
        if indices == [5, 2]:
            n = 16
        if indices == [5, 3]:
            n = 10
        if indices == [6, 0]:
            n = 12
        if indices == [6, 1]:
            n = 14
        if indices == [6, 2]:
            n = 12
        if indices == [7, 0]:
            n = 15
        if indices == [7, 1]:
            n = 13
        else:
            n = 6
        
        # EXTEND INTERVAL ON WHICH FIT IS PERFORMED TO AVOID RAPID OSCILLATIONS OF THE POLYNOMIAL AT THE INTERVAL EDGES.
        Xs_extended = []
        Xs_extended.append(Xs[0]-step)
        for X_item in Xs:
            Xs_extended.append(X_item)
        Xs_extended.append(Xs[-1]+step)
        
        Ys_extended = []
        Ys_extended.append(Ys_medians[0])
        for Y_item in Ys_medians:
            Ys_extended.append(Y_item)
        Ys_extended.append(Ys_medians[-1])
        
        # POLYNOMIAL FIT OF MEDIANS
        c, stats = P.polyfit(Xs_extended, Ys_extended, n, full=True) 
        
        x_fit = numpy.linspace(Qrod_min, Qrod_max, num=(Qrod_max-Qrod_min)/0.001, endpoint = True)
        y_fit = numpy.zeros(shape=x_fit.shape)
        for exponent in range(0, n+1):
            y_fit = y_fit + c[exponent]*(x_fit**exponent)
            
        for i in range(0, len(x_fit)):
            spot_model = [int(hRod), int(kRod), x_fit[i], y_fit[i]] # h k qRod I
            lattice_model.append(spot_model)
               
        # PLOT BRAGG ROD
        matplotlib.pyplot.scatter(cleanedList_Qrod, cleanedList_Irod, marker='o', color='c', alpha = 0.15, s=10)
        matplotlib.pyplot.plot(x_fit, y_fit, '.b-')
        myAxis = matplotlib.pyplot.gca()
        matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
        matplotlib.pyplot.axhline(y=10, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
        myAxis.set_xlim([-0.45,+0.45])
        scale = 1.1*max(cleanedList_Irod)
        myAxis.set_ylim([-0.1*scale,1*scale])
        myAxis.set_xlabel("q$_z$ (A$^{-1}$)", fontsize = 12, rotation = 'horizontal')
        matplotlib.pyplot.savefig('%s/polyFit_mergedRod_%d_%d_fit_%s.png'%(outputFolder, hRod, kRod, n))
        matplotlib.pyplot.close()
        
        # COLLECT I(qRod = 0)
        IsAtZero[rodNumber, 0]= hRod
        IsAtZero[rodNumber, 1]= kRod
        IsAtZero[rodNumber, 2]= c[0]
        
        # GENERATE braggRod OBJECT
        braggRodObject = braggRodClass.braggRod(hRod, kRod, Qrod_min, Qrod_max)
        braggRodObject.setExperimentalPoints(cleanedList_Qrod, cleanedList_Irod)
        braggRodObject.setModelPoints(x_fit, y_fit)
        braggRodObject.setModelCoefficients(c)
        joblib.dump(braggRodObject, '%s/braggRodObjects/braggRodObject_%d_%d.jbl'%(inputFolder, hRod, kRod))
                
        rodNumber = rodNumber + 1
        
    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
            
    joblib.dump(IsAtZero, '%s/intensitiesAtZero/intensitiesAtZero.jbl'%inputFolder)
    joblib.dump(lattice_model, '%s/model/lattice_model.jbl'%inputFolder)

if __name__ == "__main__":
    print "\n**** CALLING merging ****"
    mergingFunction(sys.argv[1:])    