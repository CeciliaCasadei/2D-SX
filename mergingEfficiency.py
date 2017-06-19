# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot

import makeOrbits

def mergingEfficiencyFunction(myArguments):
    
    # DEFAULTS
    outputFolder = ''
    resolutionLimit = 7.0
    
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

    print 'RUN NUMBER', runNumber
    if outputFolder == '':
        outputFolder = './Output_r%s/transformAndScale'%runNumber        
    print 'OUTPUT FOLDER ', outputFolder 
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
    myList_scaled   = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    myList_unscaled = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    
    # DEFINE ORBIT INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    orbitIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        orbitIndices.append(orbit_label)        
    print '%d Orbits'%len(orbitIndices)
    
    cellSize = 62.45
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I       

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
        
        avg_Intensity_unscaled = numpy.average(I_unscaled_cleaned)
        avg_Intensity_scaled   = numpy.average(I_scaled_cleaned)
        std_Intensity_unscaled = numpy.std(I_unscaled_cleaned)
        std_Intensity_scaled   = numpy.std(I_scaled_cleaned)
        ratio = (avg_Intensity_unscaled/std_Intensity_unscaled)  /  (avg_Intensity_scaled/std_Intensity_scaled)
        
        ratios.append(ratio)
        qs.append(q_2D)
    
    qs_plot = []
    ratios_plot = []   
    x_errors = []
    y_errors = []
    n_bins = 11
    bins = numpy.linspace(min(qs), max(qs), n_bins)
    print bins
    for i in range(0, len(bins)-1):
        left_q = bins[i]
        right_q = bins[i+1]
        x_error = right_q - left_q
        #q_bin     = [qs[i]     for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
        ratio_bin = [ratios[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
        q_mid_bin = (right_q+left_q)/2
        ratio_avg_bin = numpy.average(ratio_bin)
        qs_plot.append(q_mid_bin)
        ratios_plot.append(ratio_avg_bin)
        x_errors.append(x_error/2)
        y_errors.append(numpy.std(ratio_bin))
    print x_errors
    print y_errors
        
    matplotlib.pyplot.figure()
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=18, length=5, pad=7)    
    matplotlib.pyplot.scatter(qs, ratios, s=3)
    matplotlib.pyplot.scatter(qs_plot, ratios_plot, s=130, facecolors='m', edgecolors='m', alpha = 0.25)
    matplotlib.pyplot.errorbar(qs_plot, ratios_plot, xerr=x_errors, yerr=y_errors, capsize=0, ls='none', color='m', elinewidth=1)
    matplotlib.pyplot.gca().axhline(y=1, color='k')
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 24, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$\frac{ \delta(I^{\rm resc}) / \left < I^{\rm resc} \right > }{ \delta(I) / \left < I \right > } $", fontsize = 24, rotation = 'vertical', labelpad=20)
    matplotlib.pyplot.gca().set_ylim([0, +1.1*max(ratios)])    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/scaling_efficiency_%d_bins_IoverSig_ratios.png'%(outputFolder, n_bins))
    matplotlib.pyplot.savefig('%s/scaling_efficiency_%d_bins_IoverSig_ratios.pdf'%(outputFolder, n_bins))
        
        
                     
if __name__ == "__main__":
    print "\n**** CALLING mergingEfficiency ****"
    mergingEfficiencyFunction(sys.argv[1:])    