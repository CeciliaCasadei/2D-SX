# -*- coding: utf-8 -*-
import numpy
import scipy.optimize
import scipy.interpolate

import matplotlib.pyplot

import simulate_resolution


def prepareFile():
    directory = './F_calc_phenix'
    inputFile  = ('%s/f_model_edited.txt'%directory)
    outputFile = ('%s/f_model.txt'%directory)
    inputFile_read = open(inputFile, 'r')
    outputFile_write = open(outputFile, 'w')
    inputFile_read = list(inputFile_read)
    for i in range(0, len(inputFile_read)):
        if i%5 == 0:
            line = inputFile_read[i]
            splitLine = line.split()
            print splitLine
            h     = int(splitLine[0])
            k     = int(splitLine[1])
            l     = int(splitLine[2])
            Fobs  = float(splitLine[3])
            Fcalc = float(splitLine[8])
            outputFile_write.write('%6d%6d%6d%10.2f%10.2f\n'%(h,
                                                              k,
                                                              l,
                                                              Fobs,
                                                              Fcalc))
    
    outputFile_write.close()
    
    
def twoD_Gaussian_simple((x, y), amplitude, sigma_x, sigma_y):  
    xo = 0
    yo = 0
    g = amplitude * numpy.exp( - ((x-xo)**2)/(2*sigma_x**2) 
                               - ((y-yo)**2)/(2*sigma_y**2)   )
    return g.ravel()
    
def anisotropy():
    directory = './F_calc_phenix'
    cellSize  = 62.45     # A
    cAxis     = 180       # A
    cStar     = (2*numpy.pi)/cAxis
    
    qRods  = []
    q2Ds   = []
    Fobss  = []
    Fcalcs = []
    
    inputFile = ('%s/f_model.txt'%directory)
    inputFile_read = open(inputFile, 'r')
    for line in inputFile_read:
        splitLine = line.split()
        h     = int(splitLine[0])
        k     = int(splitLine[1])
        l     = int(splitLine[2])
        Fobs  = float(splitLine[3])
        Fcalc = float(splitLine[4])
        
        qRod = l*cStar
        d_2D = simulate_resolution.resolution(cellSize, h, k, 0.0)
        q_2D = (2*numpy.pi)/d_2D
        
        qRods.append(qRod)
        q2Ds.append(q_2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
     
    # For easier gaussian fitting.
    for i in range(0, len(Fobss)):
        qRod = qRods[i]
        q2D = - q2Ds[i]
        Fobs = Fobss[i]
        Fcalc = Fcalcs[i]
        
        qRods.append(qRod)
        q2Ds.append(q2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
        
    Fobss   = numpy.asarray(Fobss)
    Fcalcs  = numpy.asarray(Fcalcs)
    qRods   = numpy.asarray(qRods)
    q2Ds    = numpy.asarray(q2Ds)
    
    # PLOT sigma_x, sigma_y as a function of 2_2D, q_z.
    matplotlib.pyplot.figure()  
    matplotlib.pyplot.title(r'$F_{\rm obs}$', fontsize = 22)      
    matplotlib.pyplot.scatter(q2Ds, qRods, c=Fobss, 
                              #vmin=1, vmax=4, 
                              s=25, edgecolors='none')
    axes = matplotlib.pyplot.gca()
    matplotlib.pyplot.cool()
    matplotlib.pyplot.colorbar()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/Fobs.png'%(directory))
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.title(r'$F_{\rm calc}$', fontsize = 22)          
    matplotlib.pyplot.scatter(q2Ds, qRods, c=Fcalcs, 
                              #vmin=1, vmax=4, 
                              s=25, edgecolors='none')
    axes = matplotlib.pyplot.gca()
    matplotlib.pyplot.cool()
    matplotlib.pyplot.colorbar()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/Fcalc.png'%(directory))
    matplotlib.pyplot.close()
    
    #MODEL 
    initial_guess = (max(Fobss), 1.0, 1.0)
    popt_Fobss, \
    pcov_Fobss = scipy.optimize.curve_fit(twoD_Gaussian_simple, 
                                         (q2Ds, qRods), 
                                          Fobss, 
                                          p0=initial_guess)
                                          
    xx, yy = numpy.meshgrid(numpy.linspace(-1.2, 1.2, 1001),
                            numpy.linspace(-1.2, 1.2, 1001))

    zz_Fobs = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_Fobs[i,j] = twoD_Gaussian_simple((xx[i,j], yy[i,j]), 
                                                 popt_Fobss[0], 
                                                 popt_Fobss[1], 
                                                 popt_Fobss[2])  
    print popt_Fobss
    print popt_Fobss[1]/popt_Fobss[2]
    
    # PLOT THE MODELING FUNCTION
    fig, axes = matplotlib.pyplot.subplots(figsize=(5,5), facecolor='c')
    matplotlib.pyplot.title(r'F$_{\rm obs}$')
    matplotlib.pyplot.pcolor(xx, yy, zz_Fobs)
    matplotlib.pyplot.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    matplotlib.pyplot.colorbar()
    #axes = matplotlib.pyplot.gca()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=1, edgecolors='none')
    #matplotlib.pyplot.axis('equal')
    axes.set_aspect('equal')
    fig.canvas.draw()
    matplotlib.pyplot.tight_layout()    
    matplotlib.pyplot.savefig('%s/Fit_Fobs.png'%(directory), 
                              dpi=96*4)
    matplotlib.pyplot.close()
    
    #MODEL 
    popt_Fcalcs, \
    pcov_Fcalcs = scipy.optimize.curve_fit(twoD_Gaussian_simple, 
                                          (q2Ds, qRods), 
                                           Fcalcs, 
                                           p0=initial_guess)
   
    zz_Fcalc = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_Fcalc[i,j] = twoD_Gaussian_simple((xx[i,j], yy[i,j]), 
                                                  popt_Fcalcs[0], 
                                                  popt_Fcalcs[1], 
                                                  popt_Fcalcs[2])  
                                                          
    print popt_Fcalcs
    print popt_Fcalcs[1]/popt_Fcalcs[2]
    
    # PLOT THE MODELING FUNCTION
    fig, axes = matplotlib.pyplot.subplots(figsize=(5,5), facecolor='c')
    matplotlib.pyplot.title(r'F$_{\rm calc}$')
    matplotlib.pyplot.pcolor(xx, yy, zz_Fcalc)
    matplotlib.pyplot.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    matplotlib.pyplot.colorbar()
    #axes = matplotlib.pyplot.gca()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=1, edgecolors='none') 
    #matplotlib.pyplot.axis('equal')
    axes.set_aspect('equal')
    fig.canvas.draw()
    matplotlib.pyplot.tight_layout()    
    matplotlib.pyplot.savefig('%s/Fit_Fcalc.png'%(directory), 
                              dpi=96*4)
    matplotlib.pyplot.close()   


    
if __name__ == "__main__":
    print "\n**** CALLING anisotropyCheck ****"    
    anisotropy()
    
