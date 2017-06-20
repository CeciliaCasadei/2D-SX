# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

import imageSums_utilities

runNumber = '0127'
N_lattices_list = [586, 100, 10]
N_sampling_list = [10,  10,  10]
bin_size = 22


# RECIPROCAL CELL
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I                                                  # A^(-1)

matplotlib.pyplot.figure(figsize=(9, 11))   


colors = ['b', 'm', 'c']
for index in range(0, len(N_lattices_list)):
    print 'index', index
    N_lattices = N_lattices_list[index]
    N_sampling = N_sampling_list[index]
    
    ### CALCULATE CC_half OVERALL AND IN RESOLUTION BINS    
    CCs_samplings_overall = []      # N samplings
    
    CCs_samplings_bins = []         # N samplings x N bins
    qs_samplings_bins = []
    
    
    for sampling in range(0, N_sampling):
        print sampling
        
        I1_sampling = []
        I2_sampling = []
        
        CC_sampling_overall = 0
        
        CC_sampling_bins = []
        q_sampling_bins = []
            
        fOpen = open('./Output_r%s/Output_imageSums_moduleDisplacements_CChalf/imageSums_displaced_modules_%d_lattices_halves_%d.txt'%(runNumber, N_lattices, sampling))
        # h k I1_ellipse_fixed_sigmas I2_ellipse_fixed_sigmas N1 N2
        # Orbits from low to high resolution
        
        n_line = 0
        I1_bin = []
        I2_bin = []
        q_bin = []
        
        for l in fOpen:                
            if l[0] != '#':
                n_line = n_line + 1
                l_split = l.split()
                
                h = int(l_split[0])
                k = int(l_split[1])
                
                reciprocalVector = [h, k]*reciprocalCellRows
                q_x = reciprocalVector[0,0]         # A^(-1)
                q_y = reciprocalVector[0,1]         # A^(-1)
                q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
                
                I1 = float(l_split[2])
                I2 = float(l_split[3])
                I1_sampling.append(I1)
                I2_sampling.append(I2)
                
                I1_bin.append(I1)
                I2_bin.append(I2)
                q_bin.append(q)
                
            if l[0] != '#' and (n_line % bin_size) == 0:
                CC_bin = imageSums_utilities.Correlate(I1_bin, I2_bin)
                CC_sampling_bins.append(CC_bin)
                
                q_bin_avg = numpy.average(q_bin)
                q_sampling_bins.append(q_bin_avg)
                
                I1_bin = []
                I2_bin = [] 
                q_bin = []
                
        fOpen.close()
        
        CC_sampling_overall = imageSums_utilities.Correlate(I1_sampling, I2_sampling)
        CCs_samplings_overall.append(CC_sampling_overall)
        
        CCs_samplings_bins.append(CC_sampling_bins)
        qs_samplings_bins.append(q_sampling_bins)
        
        
    
    CC_overall_avg = numpy.average(CCs_samplings_overall)
    print "N_LATTICES: %6d   N_SAMPLINGS: %6d   OVERALL_CC_HALF: %.4f"%(N_lattices, N_sampling, CC_overall_avg)
    
    CCs_samplings_bins = numpy.matrix(CCs_samplings_bins)
    
    qs_samplings_bins = numpy.matrix(qs_samplings_bins)
    
    bin_avg_CCs = []
    bin_avg_qs  = []
    for j in range(0, CCs_samplings_bins.shape[1]):     ### LOOP ON BINS
        bin_avg_CC = 0
        bin_avg_q = 0
        for i in range(0, CCs_samplings_bins.shape[0]): ### LOOP ON SAMPLINGS
            CC = CCs_samplings_bins[i, j]
            q  = qs_samplings_bins[i, j]
            bin_avg_CC = bin_avg_CC + CC
            bin_avg_q  = bin_avg_q  + q
        bin_avg_CC = bin_avg_CC / CCs_samplings_bins.shape[0]
        bin_avg_q  = bin_avg_q  / CCs_samplings_bins.shape[0]
        bin_avg_CCs.append(bin_avg_CC)
        bin_avg_qs.append(bin_avg_q)
        
    matplotlib.pyplot.plot(bin_avg_qs, bin_avg_CCs, '-o', color=colors[index], mew=0)
    
matplotlib.pyplot.gca().set_ylim([0, 1.1])
matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, length=6, pad=6) 
matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 27, rotation = 'horizontal')
matplotlib.pyplot.gca().set_ylabel(r"CC$_{1/2}$",       fontsize = 25, rotation = 'vertical')
matplotlib.pyplot.savefig('./Output_r%s/Output_imageSums_moduleDisplacements_CChalf/CChalf_imageSums_displaced_modules_%d_%d_%d_lattices_%d_%d_%d_samplings_avg.png'
                          %(runNumber, N_lattices_list[0], N_lattices_list[1], N_lattices_list[2], N_sampling_list[0], N_sampling_list[1], N_sampling_list[2]), dpi=4*96)
matplotlib.pyplot.savefig('./Output_r%s/Output_imageSums_moduleDisplacements_CChalf/CChalf_imageSums_displaced_modules_%d_%d_%d_lattices_%d_%d_%d_samplings_avg.pdf'
                          %(runNumber, N_lattices_list[0], N_lattices_list[1], N_lattices_list[2], N_sampling_list[0], N_sampling_list[1], N_sampling_list[2]), dpi=4*96)
matplotlib.pyplot.close()