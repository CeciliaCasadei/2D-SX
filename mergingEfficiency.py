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
        
        std_Intensity_unscaled = numpy.std(I_unscaled_cleaned)
        std_Intensity_scaled   = numpy.std(I_scaled_cleaned)
        ratio = std_Intensity_scaled/std_Intensity_unscaled
        
        ratios.append(ratio)
        qs.append(q_2D)
    
    qs_plot = []
    ratios_plot = []   
    bins = numpy.linspace(min(qs), max(qs), 20)
    for i in range(0, len(bins)-1):
        left_q = bins[i]
        right_q = bins[i+1]
        q_bin     = [qs[i]     for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
        ratio_bin = [ratios[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
        q_avg_bin = numpy.average(q_bin)
        ratio_avg_bin = numpy.average(ratio_bin)
        qs_plot.append(q_avg_bin)
        ratios_plot.append(ratio_avg_bin)
        
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(qs, ratios, s=3)
    matplotlib.pyplot.scatter(qs_plot, ratios_plot, s=130, facecolors='r', edgecolors='r', alpha = 0.25)
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$\frac{\sigma({\rm I}_{\rm scaled})}{\sigma({\rm I}_{\rm unscaled})} $", fontsize = 18, rotation = 'horizontal', labelpad=30)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('./scaling_efficiency_20bins.png')
        
        
                     
if __name__ == "__main__":
    print "\n**** CALLING mergingEfficiency ****"
    mergingEfficiencyFunction(sys.argv[1:])    