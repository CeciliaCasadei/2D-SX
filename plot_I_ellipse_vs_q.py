# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

runNumber = '0127'

outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I                     
        
data = open('%s/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse_x0_y0.txt'%outputFolder, 'r')
Is_ellipse = []
qs = []
for data_line in data:
    h = int(data_line.split()[0])
    k = int(data_line.split()[1])
    I_ellipse = float(data_line.split()[4])
    
    reciprocalVector = [h, k]*reciprocalCellRows
    q_x = reciprocalVector[0,0]         # A^(-1)
    q_y = reciprocalVector[0,1]         # A^(-1)
    q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
    
    qs.append(q)
    Is_ellipse.append(I_ellipse)
    
    
# BINNING    
qs_bins = []
Is_ellipse_bins = [] 


bins = numpy.linspace(min(qs), max(qs), 25)
for i in range(0, len(bins)-1):
    left_q  = bins[i]
    right_q = bins[i+1]
    q_bin               = [qs[i]         for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]    
    I_ellipse_bin       = [Is_ellipse[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
   
    print len(q_bin)
    #q_avg_bin = numpy.average(q_bin)
    q_mid_bin = (left_q+right_q)/2
    I_ellipse_avg_bin = numpy.average(I_ellipse_bin)
    
    qs_bins.append(q_mid_bin)
    Is_ellipse_bins.append(I_ellipse_avg_bin)
    

print len(qs), len(Is_ellipse)   
matplotlib.pyplot.scatter(qs, Is_ellipse, s=3)
matplotlib.pyplot.scatter(qs_bins, Is_ellipse_bins, s=100, facecolors='c', edgecolors='none', alpha = 0.40)
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_ylim([0.05, 500])
matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=12, length=5, pad=7)    
matplotlib.pyplot.gca().set_xlabel(r'$q$ ($\AA^{-1}$)', fontsize = 18)
matplotlib.pyplot.gca().set_ylabel(r"I (n photons)", rotation = 'vertical', fontsize = 18)
matplotlib.pyplot.tight_layout()
matplotlib.pylab.savefig('%s/I_ellipse_fixed_sigmas_vs_q.png'%outputFolder, dpi=96*4)
matplotlib.pylab.savefig('%s/I_ellipse_fixed_sigmas_vs_q.pdf'%outputFolder, dpi=96*4)