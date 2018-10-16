# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import numpy



import braggRodClass
import makeOrbits


def mergingFunction(myArguments):
    input_str = '--inputFolder <inputFolder> --resolutionLimit <resolutionLimit>'
    
    # DEFAULTS
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

    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
        
    print '%d Rods'%len(rodIndices)   
                  
    # PLOT I vs QROD, FOR EVERY ROD
    for indices in rodIndices:
        print '\n\n************************\nRod: %s'%indices
        hRod = indices[0]
        kRod = indices[1]
        Qrod_vector  = []
        Irod_vector  = []
        dIrod_vector = []
        
        # FOR EVERY ROD, COLLECT (QROD, I, dI) POINTS FROM ALL RUNS 
        for runNumber in runNumbers:
            print 'Extracting (qRod, I) points from run %s'%runNumber
            listDirectory = '%s/spotsMatricesList-Scaled-r%s'%(inputFolder, 
                                                               runNumber)
            myList = joblib.load('%s/r%s_scaledSpotsMatricesList_error.jbl'
                                  %(listDirectory, runNumber))
        
            for latticeMatrix in myList: 
                latticeMatrix = numpy.asarray(latticeMatrix) # h 
                                                             # k 
                                                             # qRod 
                                                             # I (Photons, LP-corrected, scaled)
                                                             # flag
                                                             # i (unassembled matrix)
                                                             # j (unassembled matrix)
                                                             # dI
                if latticeMatrix[0, 4] == 1:
                    for spot in latticeMatrix:
                        h = spot[0]
                        k = spot[1]
                        if ((h == hRod       and k == kRod) or 
                            (h == -hRod-kRod and k == hRod) or 
                            (h == kRod       and k == -hRod-kRod)):
                            Qrod_vector.append(spot[2])
                            Irod_vector.append(spot[3])
                            dIrod_vector.append(spot[7])
                            
                            
                        if ((h == -hRod     and k == -kRod) or 
                            (h == hRod+kRod and k == -hRod) or 
                            (h == -kRod     and k == hRod+kRod)):                            
                            Qrod_vector.append(-spot[2])
                            Irod_vector.append(spot[3])
                            dIrod_vector.append(spot[7])
                        
                        #print spot[3], spot[7]
                            
         
        # REMOVE NAN VALUES (Coming from Is whose integration area is not in one detector module)
        cleanedList_Irod  = [Irod_vector[i]  for i in range(0, len(Irod_vector)) 
                                             if not numpy.isnan(Irod_vector[i])]
        cleanedList_Qrod  = [Qrod_vector[i]  for i in range(0, len(Irod_vector)) 
                                             if not numpy.isnan(Irod_vector[i])]
        cleanedList_dIrod = [dIrod_vector[i] for i in range(0, len(Irod_vector)) 
                                             if not numpy.isnan(Irod_vector[i])]
                
        Qrod_min = min(cleanedList_Qrod)
        Qrod_max = max(cleanedList_Qrod)


               
        # PLOT BRAGG ROD
        matplotlib.pyplot.scatter(cleanedList_Qrod, 
                                  cleanedList_Irod, 
                                  marker='o', 
                                  color='c', 
                                  alpha = 0.15, 
                                  s=10)
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
    

    
        # GENERATE braggRod OBJECT
        braggRodObject = braggRodClass.braggRod(hRod, 
                                                kRod, 
                                                Qrod_min, 
                                                Qrod_max)
        braggRodObject.setExperimentalPoints(cleanedList_Qrod, 
                                             cleanedList_Irod)
        
        braggRodObject.setExperimentalErrors(cleanedList_dIrod)

        joblib.dump(braggRodObject, '%s/braggRodObject_%d_%d.jbl'
                                     %(rodsFolder, hRod, kRod))
                                     
   

if __name__ == "__main__":
    print "\n**** CALLING merging ****"
    mergingFunction(sys.argv[1:])    