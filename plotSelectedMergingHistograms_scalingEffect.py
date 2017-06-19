# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot


def plotScalingHistograms(myArguments):
    
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
    
    
    
    # DEFINE ORBIT INDICES 
    orbitIndices = [[1, 1], [3, 4]]        
    
    cellSize = 62.45
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    
    # 4 SUBPLOTS    
    myFigure, a = matplotlib.pyplot.subplots(2, 2)    
    a = a.ravel()
    
    myLabels = ['a', 'b', 'c', 'd']
    
    # FOR EACH ORBIT, COLLECT P3-RELATED POINTS   
    for idx, ax in enumerate(a):   
        indices = orbitIndices[idx/2]
           
        hRod = indices[0]
        kRod = indices[1]
        reciprocalVector = [hRod, kRod]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
        resolution = 2* numpy.pi / q_2D     # A
        
        
        print idx, indices
        
        if idx%2 == 0:
            print 'unscaled'
            myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
            
        else:
            print 'scaled'
            
            # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
            myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
        
        Iorbit_vector = []
            
            
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
        ax.axvline(x=avg_Intensity,                  ymin = 0, ymax = 1, color='r')
        ax.axvline(x=avg_Intensity-stddev_Intensity, ymin = 0, ymax = 1, color='r', linestyle=':')
        ax.axvline(x=avg_Intensity+stddev_Intensity, ymin = 0, ymax = 1, color='r', linestyle=':')
        
        ax.tick_params(axis='both', which='major', labelsize=10, length=5)          
        ax.text(0.96, 0.94, "{(%d, %d)}\n%.1f $\AA$\nN$_{\mathrm{tot}}$ = %d\n(%.1f $\pm$ %.1f) photons"%(hRod, kRod, resolution, len(Iorbit_vector_cleaned), avg_Intensity, stddev_Intensity), horizontalalignment='right', verticalalignment='top', fontsize=10, transform = ax.transAxes) 
        ax.text(0.03, 0.03, "(%s)"%myLabels[idx], horizontalalignment='left', verticalalignment='bottom', fontsize=16, transform = ax.transAxes) 
        ax.set_xlabel(r"$I$ (n photons)", fontsize = 14, rotation = 'horizontal')
        ax.set_ylabel(r"$n_{\rm obs}$", fontsize = 14, rotation = 'horizontal', labelpad = 14)
        ax.set_xlim([-1*avg_Intensity, +5*avg_Intensity]) # was 3, 5
        ax.set_ylim([0, +1.2*max(n)])

    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/Histograms/4xmergingHistograms_scalingEffect_r%s.png'%(outputFolder, runNumber), dpi = 4*96, facecolor = 'w')
    matplotlib.pyplot.savefig('%s/Histograms/4xmergingHistograms_scalingEffect_r%s.pdf'%(outputFolder, runNumber), facecolor = 'w')
    matplotlib.pyplot.close()
        
                     
if __name__ == "__main__":
    print "\n**** CALLING plotSelectedMergingHistograms_cmp_scalingEffects ****"
    plotScalingHistograms(sys.argv[1:])    