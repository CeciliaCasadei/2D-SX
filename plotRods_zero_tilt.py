# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot


import makeOrbits

def plotRodsFunction(myArguments):
    
    # DEFAULTS
    outputFolder = ''
    resolutionLimit = 4.0
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "outputFolder="])
    except getopt.GetoptError:
        print 'Usage: python plotRods.py --runNumber <runNumber> --outputFolder <outputFolder>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRods.py --runNumber <runNumber> --outputFolder <outputFolder>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--outputFolder":
            outputFolder = value

    print runNumber
    if outputFolder == '':
        outputFolder = './Output_r%s/transformAndScale'%runNumber
        
    print 'OUTPUT FOLDER ', outputFolder
    outputFile = open('%s/mergedIntensities_h_k_qRod_l_avgI_sigmaI.txt'%outputFolder, 'w')
    
    # FOLDERS
    if not os.path.exists('%s/Rods'%outputFolder):
        os.mkdir('%s/Rods'%outputFolder)            
    if not os.path.exists('%s/Histograms'%outputFolder):
        os.mkdir('%s/Histograms'%outputFolder)
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
    myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    
    # DEFINE ROD INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
        
    print '%d Rods'%len(rodIndices)
    
#    cellSize = 62.45
#    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
#    reciprocalCellRows = 2 * numpy.pi * directCell.I    
#    
    #FOR EACH ROD, COLLECT P3-RELATED POINTS AND FRIEDEL MATES              
    for indices in rodIndices:
        hRod = indices[0]
        kRod = indices[1]
        Irod_vector = []
        Qrod_vector = []
        Irod_vector_I = []
        Irod_vector_i = []
        Qrod_vector_I = []
        Qrod_vector_i = []
        
     #   reciprocalVector = [hRod, kRod]*reciprocalCellRows
      #  q_x = reciprocalVector[0,0]         # A^(-1)
       # q_y = reciprocalVector[0,1]         # A^(-1)
        #q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        
        for latticeMatrix in myList:   
            latticeMatrix = numpy.asarray(latticeMatrix)
            if latticeMatrix[0, 4] == 0:                                         # Check lattice flag
                continue
            for spot in latticeMatrix:
                h_transformed = spot[0]
                k_transformed = spot[1]
                if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                    Irod_vector.append(spot[3])
                    Qrod_vector.append(spot[2])
                    
                    Irod_vector_I.append(spot[3])
                    Qrod_vector_I.append(spot[2])
                    
                if (h_transformed == -hRod and k_transformed == -kRod) or (h_transformed == hRod+kRod and k_transformed == -hRod) or (h_transformed == -kRod and k_transformed == hRod+kRod):
                    Irod_vector.append(spot[3])
                    Qrod_vector.append(-spot[2])  ### NB minus!!!
                    
                    Irod_vector_i.append(spot[3])
                    Qrod_vector_i.append(spot[2]) ### NB plus!!!

        # PLOT ROD
        matplotlib.pyplot.scatter(Qrod_vector, Irod_vector, marker='.', color='b', alpha=0.15)
        if not os.path.exists('%s/Rods'%outputFolder):
            os.mkdir('%s/Rods'%outputFolder)
        matplotlib.pyplot.savefig('%s/Rods/r%s_rod_%d_%d.png'%(outputFolder, runNumber, hRod, kRod))
        matplotlib.pyplot.close()
        
        # LOG RESULTS
        Irod_vector_I_cleaned = [Irod_vector_I[i] for i in range(0, len(Irod_vector_I)) if not numpy.isnan(Irod_vector_I[i])]
        Irod_vector_i_cleaned = [Irod_vector_i[i] for i in range(0, len(Irod_vector_i)) if not numpy.isnan(Irod_vector_i[i])]
        Qrod_vector_I_cleaned = [Qrod_vector_I[i] for i in range(0, len(Qrod_vector_I)) if not numpy.isnan(Irod_vector_I[i])]
        Qrod_vector_i_cleaned = [Qrod_vector_i[i] for i in range(0, len(Qrod_vector_i)) if not numpy.isnan(Irod_vector_i[i])]
        
        avg_Intensity_I = numpy.average(Irod_vector_I_cleaned)
        avg_Intensity_i = numpy.average(Irod_vector_i_cleaned)
        
        avg_Qrod_I = numpy.average(Qrod_vector_I_cleaned)
        avg_Qrod_i = numpy.average(Qrod_vector_i_cleaned)
        
        stddev_Intensity_I = numpy.std(Irod_vector_I_cleaned)
        stddev_Intensity_i = numpy.std(Irod_vector_i_cleaned)
        
        outputFile.write('\n%3d %3d %8.4f %3d %7.2f %7.2f'%( hRod,  kRod, avg_Qrod_I, int( avg_Qrod_I*504.5/(2*numpy.pi)), avg_Intensity_I, stddev_Intensity_I))
        outputFile.write('\n%3d %3d %8.4f %3d %7.2f %7.2f'%(-hRod, -kRod, avg_Qrod_i, int( avg_Qrod_i*504.5/(2*numpy.pi)), avg_Intensity_i, stddev_Intensity_i))
                
        # PLOT ORBIT HISTOGRAMS
        Irod_vector_I_cleaned = numpy.ravel(Irod_vector_I_cleaned)
        Irod_vector_i_cleaned = numpy.ravel(Irod_vector_i_cleaned)
        
        matplotlib.pyplot.figure(facecolor = 'w')
        matplotlib.pyplot.hist(Irod_vector_I_cleaned)
        matplotlib.pyplot.xlabel(r"I$_{\rm scaled}$ (photons)")
        matplotlib.pyplot.ylabel(r"N$_{\rm measured}$")        
        matplotlib.pyplot.savefig('%s/Histograms/Histogram_r%s_spot_%d_%d.png'%(outputFolder, runNumber, hRod, kRod), dpi = 300, facecolor = 'w')
        matplotlib.pyplot.close()
               
        matplotlib.pyplot.figure(facecolor = 'w')
        matplotlib.pyplot.hist(Irod_vector_i_cleaned)
        matplotlib.pyplot.xlabel(r"I$_{\rm scaled}$ (photons)")
        matplotlib.pyplot.ylabel(r"N$_{\rm measured}$")
        matplotlib.pyplot.savefig('%s/Histograms/Histogram_r%s_spot_%d_%d.png'%(outputFolder, runNumber, -hRod, -kRod), dpi = 300, facecolor = 'w')
        matplotlib.pyplot.close()
        
            
    outputFile.close()        

if __name__ == "__main__":
    print "\n**** CALLING plotRods ****"
    plotRodsFunction(sys.argv[1:])    