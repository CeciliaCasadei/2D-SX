# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot


def plotOrbitsFunction(myArguments):
    
    # DEFAULTS
    outputFolder = ''
    
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
    
    # FOLDERS        
    if not os.path.exists('%s/Histograms'%outputFolder):
        os.mkdir('%s/Histograms'%outputFolder)
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
    myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    
    # DEFINE ORBIT INDICES 
    orbitIndices = [[1, 1], [7, 1], [3, 4], [2, 11]]        
    
    cellSize = 62.45
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    
    # 4 SUBPLOTS    
    myFigure, a = matplotlib.pyplot.subplots(2, 2)    
    a = a.ravel()
    
    # FOR EACH ORBIT, COLLECT P3-RELATED POINTS   
    for idx, ax in enumerate(a):        
        indices = orbitIndices[idx]
        
        hRod = indices[0]
        kRod = indices[1]
        Iorbit_vector = []
        
        reciprocalVector = [hRod, kRod]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        resolution = 2* numpy.pi / q_2D     # A
        
        for latticeMatrix in myList:   
            latticeMatrix = numpy.asarray(latticeMatrix)
            if latticeMatrix[0, 4] == 0:                                         # Check lattice flag
                continue
            for spot in latticeMatrix:
                h_transformed = spot[0]
                k_transformed = spot[1]
                if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                    Iorbit_vector.append(spot[3])
                    
        Iorbit_vector_cleaned = [Iorbit_vector[i] for i in range(0, len(Iorbit_vector)) if not numpy.isnan(Iorbit_vector[i])]       
        avg_Intensity = numpy.average(Iorbit_vector_cleaned)
        stddev_Intensity = numpy.std(Iorbit_vector_cleaned)       
        Iorbit_vector_cleaned = numpy.ravel(Iorbit_vector_cleaned)   
    
        # PLOT HISTOGRAM
        maxI = max(Iorbit_vector_cleaned)
        minI = min(Iorbit_vector_cleaned)
        step = avg_Intensity/4
        myBins = numpy.arange(minI-step, maxI+step, step)        
        (n, bins, patches) = ax.hist(Iorbit_vector_cleaned, bins=myBins)
        ax.axvline(x=0, color='k', linestyle=':')
        ax.axvline(x=avg_Intensity,                    ymin = 0, ymax = 0.9, color='r')
        ax.axvline(x=avg_Intensity-stddev_Intensity, ymin = 0, ymax = 0.9, color='r', linestyle=':')
        ax.axvline(x=avg_Intensity+stddev_Intensity, ymin = 0, ymax = 0.9, color='r', linestyle=':')
        ax.set_xlim([-6*avg_Intensity, +8*avg_Intensity]) # was 3, 5
        ax.set_ylim([0, +1.2*max(n)])
        ax.tick_params(axis='both', which='major', labelsize=8)       
        ax.text(0.03, 0.97, "{(%d, %d)}\n%.1f $\AA$"%(hRod, kRod, resolution), horizontalalignment='left', verticalalignment='top', fontsize=10, transform = ax.transAxes) 
        ax.text(0.97, 0.97, "(%.1f $\pm$ %.1f) photons"%(avg_Intensity, stddev_Intensity), horizontalalignment='right', verticalalignment='top', fontsize=10, transform = ax.transAxes) 

    matplotlib.pyplot.savefig('%s/Histograms/4xHIS_r%s.png'%(outputFolder, runNumber), dpi = 4*96, facecolor = 'w')
    matplotlib.pyplot.close()
        

                     
if __name__ == "__main__":
    print "\n**** CALLING plotSelectedMergingHistograms ****"
    plotOrbitsFunction(sys.argv[1:])    