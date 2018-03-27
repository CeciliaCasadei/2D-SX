# -*- coding: utf-8 -*-
import numpy
import scipy.optimize
import scipy.interpolate

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import simulate_resolution


def prepareFile():
    directory = './F_calc_phenix'
    inputFile  = ('%s/f_model_onlyBulkSolvent_edited.txt'%directory)
    outputFile = ('%s/f_model_onlyBulkSolvent.txt'%directory)
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
            Fcalc = float(splitLine[6])
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
    
def ratioModel((x, y), amplitude, B_x, B_y):  
    g = amplitude * numpy.exp( - (B_x*(x**2))
                               - (B_y*(y**2))  )
    return g.ravel()
    
    
def lnRatioModel((x, y), amplitude, B_x, B_y):  
    g = amplitude  - B_x*x - B_y*y
    return g.ravel()
    
    
def anisotropy():
    directory = './F_calc_phenix'
    cellSize  = 62.45     # A
    
    # Double sampling!
    #cAxis     = 180       # A
    
    # Single sampling!
    cAxis     = 90        # A
    cStar     = (2*numpy.pi)/cAxis
    
    qRods  = []
    q2Ds   = []
    Fobss  = []
    Fcalcs = []
    ratios = []
    
    inputFile = ('%s/f_model_onlyBulkSolvent.txt'%directory)
    inputFile_read = open(inputFile, 'r')
    for line in inputFile_read:
        splitLine = line.split()
        h     = int(splitLine[0])
        k     = int(splitLine[1])
        l     = int(splitLine[2])
        Fobs  = float(splitLine[3])
        Fcalc = float(splitLine[4]) # PDB MODEL AND BULK SOLVENT, NO ANISOTROPY CORRECTION
        ratio = (Fobs**2)/(Fcalc**2)
        
        qRod = l*cStar
        d2D = simulate_resolution.resolution(cellSize, h, k, 0.0)
        q2D = (2*numpy.pi)/d2D
        
        qRods.append(abs(qRod))
        q2Ds.append(q2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
        ratios.append(ratio)
           
    Fobss   = numpy.asarray(Fobss)
    Fcalcs  = numpy.asarray(Fcalcs)
    qRods   = numpy.asarray(qRods)
    q2Ds    = numpy.asarray(q2Ds)
    ratios  = numpy.asarray(ratios)
    
    #MODEL
    q2Ds_clean   = [q2Ds[i]   for i in range(0, len(q2Ds))]# if ratios[i] < 40000]
    qRods_clean  = [qRods[i]  for i in range(0, len(q2Ds))]# if ratios[i] < 40000]
    ratios_clean = [ratios[i] for i in range(0, len(q2Ds))]# if ratios[i] < 40000]
    q2Ds_clean   = numpy.asarray(q2Ds_clean)
    qRods_clean  = numpy.asarray(qRods_clean)
    ratios_clean = numpy.asarray(ratios_clean)
    nOutliers = ratios.shape[0] - ratios_clean.shape[0]
    print "N OUTLIERS: ", nOutliers

    lnRatios = numpy.log(ratios_clean)
    q2D_squared = q2Ds_clean**2
    qRod_squared = qRods_clean**2
    
    popt_lnRatios, \
    pcov_lnRatios = scipy.optimize.curve_fit(lnRatioModel, 
                                            (q2D_squared, qRod_squared), 
                                            lnRatios)
                                            
    print "ln(Ratios)"
    print popt_lnRatios
   
    xx, yy = numpy.meshgrid(numpy.linspace(0, max(q2D_squared), 101),
                            numpy.linspace(0, max(qRod_squared), 101))
                            
    zz_ratio = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_ratio[i,j] = lnRatioModel((xx[i,j], yy[i,j]), 
                                         popt_lnRatios[0], 
                                         popt_lnRatios[1], 
                                         popt_lnRatios[2]) 
    
    
    # 3D Figure
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(q2D_squared, qRod_squared, lnRatios, c='k', marker='o', s=0.5)
    mySurface = myAxes.plot_surface(xx, yy, zz_ratio, 
                                    rstride=1, cstride=1, 
                                    cmap=cm.coolwarm,
                                    linewidth=0, 
                                    antialiased=False, 
                                    alpha = 0.2)
    myFigure.colorbar(mySurface, shrink=0.5, aspect=5)    
    myAxes.zaxis.set_major_locator(LinearLocator(10))
    myAxes.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    myAxes.set_xlabel(r'q$^{2}$$_{\rm 2D}$ [A$^{-2}$]')
    myAxes.set_ylabel(r'q$^{2}$$_{\rm rod}$ [A$^{-2}$]') 
    myAxes.set_zlabel(r'I$_{\rm obs}$/I$_{\rm calc}$') 
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
    matplotlib.pyplot.tight_layout()  
    print 'Showing 3D figure.'
    #matplotlib.pyplot.show()
    matplotlib.pyplot.savefig('./F_calc_phenix/ratios.png', dpi=4*96)
    
if __name__ == "__main__":
    print "\n**** CALLING anisotropyCheck ****"  
    prepareFile()
    anisotropy()
    
