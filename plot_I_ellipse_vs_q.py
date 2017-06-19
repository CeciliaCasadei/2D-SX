# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot
outputFolder = './Output_imageSums_moduleDisplacements_sigmaFits'
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

print len(qs), len(Is_ellipse)   
matplotlib.pyplot.scatter(qs, Is_ellipse, s=3)
matplotlib.pyplot.gca().set_yscale('log')
matplotlib.pyplot.gca().set_ylim([0.01,1000])
matplotlib.pyplot.gca().set_xlabel(r'q ($\AA^{-1}$)')
matplotlib.pyplot.gca().set_ylabel(r'I$_{\rm ellipse}$ (photon counts)', rotation = 'vertical')
matplotlib.pylab.savefig('%s/I_ellipse_fixed_sigmas_vs_q.png'%outputFolder, dpi=96*4)
matplotlib.pylab.savefig('%s/I_ellipse_fixed_sigmas_vs_q.pdf'%outputFolder, dpi=96*4)