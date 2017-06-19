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
    orbitIndices = [[1, 1], [7, 1], [3, 4]]        
    
    cellSize = 62.45
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    
    # 4 SUBPLOTS    
    myFigure, a = matplotlib.pyplot.subplots(2, 2)    
    a = a.ravel()
    
    # FOR EACH ORBIT, COLLECT P3-RELATED POINTS   
    for idx, ax in enumerate(a):   
        try:
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
            ax.set_xlim([-3*avg_Intensity, +5*avg_Intensity]) # was 3, 5
            ax.set_ylim([0, +1.2*max(n)])
            ax.tick_params(axis='both', which='major', labelsize=8)       
            ax.text(0.03, 0.97, "{(%d, %d)}\n%.1f $\AA$\nN$_{\mathrm{tot}}$ = %d"%(hRod, kRod, resolution, len(Iorbit_vector_cleaned)), horizontalalignment='left', verticalalignment='top', fontsize=10, transform = ax.transAxes) 
            ax.text(0.97, 0.97, "(%.1f $\pm$ %.1f) photons"%(avg_Intensity, stddev_Intensity), horizontalalignment='right', verticalalignment='top', fontsize=10, transform = ax.transAxes) 
            ax.set_xlabel(r"I (photons)", fontsize = 10, rotation = 'horizontal')
            ax.set_ylabel(r"N", fontsize = 10, rotation = 'horizontal', labelpad = 10)
        except:
            print 'Plotting merging efficiency'
            
            # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
            myList_scaled   = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
            myList_unscaled = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
            
            # DEFINE ORBIT INDICES 
            orbits = makeOrbits.makeOrbitsFunction(7)
            orbitIndices = []
            for orbit in orbits:
                orbit_label = orbit.label
                orbitIndices.append(orbit_label)        
            print '%d Orbits'%len(orbitIndices)
            
            ratios = []
            qs = []
            #FOR EACH ORBIT, COLLECT P3-RELATED POINTS
            for indices in orbitIndices:
                hOrbit = indices[0]
                kOrbit = indices[1]
                print hOrbit, kOrbit
                
                I_scaled   = []
                I_unscaled = []
                
                reciprocalVector = [hOrbit, kOrbit]*reciprocalCellRows
                q_x = reciprocalVector[0,0]         # A^(-1)
                q_y = reciprocalVector[0,1]         # A^(-1)
                q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
                
                flags = []
                for latticeMatrix_scaled in myList_scaled:   
                    latticeMatrix_scaled = numpy.asarray(latticeMatrix_scaled)
                    flags.append(latticeMatrix_scaled[0, 4])
                    if latticeMatrix_scaled[0, 4] == 0:                                         # Check lattice flag
                        continue
                    for spot in latticeMatrix_scaled:
                        h_transformed = spot[0]
                        k_transformed = spot[1]
                        if (h_transformed == hOrbit and k_transformed == kOrbit) or (h_transformed == -hOrbit-kOrbit and k_transformed == hOrbit) or (h_transformed == kOrbit and k_transformed == -hOrbit-kOrbit):
                            I_scaled.append(spot[3])
                print len(flags) 
                
                flag_idx = 0           
                for latticeMatrix_unscaled in myList_unscaled:  
                    flag_idx = flag_idx + 1
                    latticeMatrix_unscaled = numpy.asarray(latticeMatrix_unscaled)
                    if flags[flag_idx-1] == 0:
                        continue
                    
                    for spot in latticeMatrix_unscaled:
                        h_transformed = spot[0]
                        k_transformed = spot[1]
                        if (h_transformed == hOrbit and k_transformed == kOrbit) or (h_transformed == -hOrbit-kOrbit and k_transformed == hOrbit) or (h_transformed == kOrbit and k_transformed == -hOrbit-kOrbit):
                            I_unscaled.append(spot[3])
                
                print len(I_scaled), len(I_unscaled)                  
                I_unscaled_cleaned = [I_unscaled[i] for i in range(0, len(I_unscaled)) if not numpy.isnan(I_unscaled[i])]
                I_scaled_cleaned   = [I_scaled[i]   for i in range(0, len(I_scaled))   if not numpy.isnan(I_scaled[i])] 
                
                std_Intensity_unscaled = numpy.std(I_unscaled_cleaned)
                std_Intensity_scaled   = numpy.std(I_scaled_cleaned)
                ratio = std_Intensity_scaled/std_Intensity_unscaled
                
                ratios.append(ratio)
                qs.append(q_2D)
            
            qs_plot = []
            ratios_plot = []   
            bins = numpy.linspace(min(qs), max(qs), 11)
            for i in range(0, len(bins)-1):
                left_q = bins[i]
                right_q = bins[i+1]
                q_bin     = [qs[i]     for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
                ratio_bin = [ratios[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
                q_avg_bin = numpy.average(q_bin)
                ratio_avg_bin = numpy.average(ratio_bin)
                qs_plot.append(q_avg_bin)
                ratios_plot.append(ratio_avg_bin)
                
            ax.scatter(qs, ratios, s=3)
            ax.scatter(qs_plot, ratios_plot, s=100, facecolors='r', edgecolors='r', alpha = 0.25)
            ax.tick_params(axis='both', which='major', labelsize=8)    
            ax.set_xlabel(r"q ($\AA^{-1}$)", fontsize = 10, rotation = 'horizontal')
            ax.set_ylabel(r"$\frac{\sigma({\rm I}_{\rm s})}{\sigma({\rm I}_{\rm u})} $", fontsize = 10, rotation = 'horizontal', labelpad = 10)
            
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/Histograms/3xmergingHistograms_mergingEfficiency_r%s.png'%(outputFolder, runNumber), dpi = 4*96, facecolor = 'w')
    matplotlib.pyplot.savefig('%s/Histograms/3xmergingHistograms_mergingEfficiency_r%s.pdf'%(outputFolder, runNumber), facecolor = 'w')
    matplotlib.pyplot.close()
        
                     
if __name__ == "__main__":
    print "\n**** CALLING plotSelectedMergingHistograms ****"
    plotOrbitsFunction(sys.argv[1:])    