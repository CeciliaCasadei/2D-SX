# -*- coding: utf-8 -*-
import numpy
import scipy.optimize
import scipy.interpolate

import matplotlib
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

import simulate_resolution


def prepareFile(directory, inFile, outFile):
    
    inputFile  = '%s/%s'%(directory, inFile)
    outputFile = '%s/%s'%(directory, outFile)
    
    inputFile_read   = open(inputFile, 'r')
    outputFile_write = open(outputFile, 'w')
    inputFile_read   = list(inputFile_read)
    
    for i in range(0, len(inputFile_read)):
        if i%5 == 0:
            line = inputFile_read[i]
            splitLine = line.split()
            
            h       = int(splitLine[0])
            k       = int(splitLine[1])
            l       = int(splitLine[2])
            Fobs    = float(splitLine[3])
            sigFobs = float(splitLine[4])
            Fcalc   = float(splitLine[6]) # Include bulk solvent
            outputFile_write.write('%6d%6d%6d%10.2f%10.2f%10.2f\n'%(h,
                                                                    k,
                                                                    l,
                                                                    Fobs,
                                                                    sigFobs,
                                                                    Fcalc))
    
    outputFile_write.close()

def prepareFile_noSigF(directory, inFile, outFile):
    
    inputFile  = '%s/%s'%(directory, inFile)
    outputFile = '%s/%s'%(directory, outFile)
    
    inputFile_read   = open(inputFile, 'r')
    outputFile_write = open(outputFile, 'w')
    inputFile_read   = list(inputFile_read)
    
    for i in range(0, len(inputFile_read)):
        if i%4 == 0:
            line = inputFile_read[i]
            splitLine = line.split()
            
            h       = int(splitLine[0])
            k       = int(splitLine[1])
            l       = int(splitLine[2])
            Fobs    = float(splitLine[3])
            Fcalc   = float(splitLine[5]) # Include bulk solvent
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
    
    
def anisotropy(directory, inFile, cellSize, cAxis):
    
    cStar     = (2*numpy.pi)/cAxis
    
    qRods    = []
    q2Ds     = []
    Fobss    = []
    sigFobss = []
    Fcalcs   = []
    ratios   = []
    
    inputFile = ('%s/%s'%(directory, inFile))
    inputFile_read = open(inputFile, 'r')
    for line in inputFile_read:
        splitLine = line.split()
        
        h       = int(splitLine[0])
        k       = int(splitLine[1])
        l       = int(splitLine[2])
        Fobs    = float(splitLine[3])
        sigFobs = float(splitLine[4])
        Fcalc   = float(splitLine[5]) # PDB MODEL AND BULK SOLVENT, NO ANISOTROPY CORRECTION
        ratio   = (Fobs**2)/(Fcalc**2)
        
        if h == 0 and k == 0:
            continue
        
        qRod = l*cStar
        d2D  = simulate_resolution.resolution(cellSize, h, k, 0.0)
        q2D  = (2*numpy.pi)/d2D
        
        qRods.append(abs(qRod))
        q2Ds.append(q2D)
        Fobss.append(Fobs)
        sigFobss.append(sigFobs)
        Fcalcs.append(Fcalc)
        ratios.append(ratio)
           
    qRods    = numpy.asarray(qRods)
    q2Ds     = numpy.asarray(q2Ds)
    Fobss    = numpy.asarray(Fobss)
    sigFobss = numpy.asarray(sigFobss)
    Fcalcs   = numpy.asarray(Fcalcs)
    ratios   = numpy.asarray(ratios)

    lnRatios     = numpy.log(ratios)
    q2D_squared  = q2Ds**2
    qRod_squared = qRods**2
    
    popt_lnRatios, \
    pcov_lnRatios = scipy.optimize.curve_fit(lnRatioModel, 
                                             (q2D_squared, qRod_squared), 
                                             lnRatios)
                                            
    lnA    = popt_lnRatios[0]
    dB_2D  = popt_lnRatios[1]
    dB_rod = popt_lnRatios[2]
   
    print 'Parameters: ', lnA, dB_2D, dB_rod
    
    xx, yy = numpy.meshgrid(numpy.linspace(0, max(q2D_squared),  101),
                            numpy.linspace(0, max(qRod_squared), 101))
                            
    zz_ratio = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_ratio[i,j] = lnRatioModel((xx[i,j], yy[i,j]), 
                                         lnA, 
                                         dB_2D, 
                                         dB_rod)
    
    # 3D Figure
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(q2D_squared, 
                   qRod_squared, 
                   lnRatios, 
                   c='k', 
                   marker='o', 
                   s=0.5)
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
    myAxes.set_zlabel(r'ln(I$_{\rm obs}$/I$_{\rm calc}$)') 
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
    matplotlib.pyplot.tight_layout()  
    print 'Showing 3D figure.'
    matplotlib.pyplot.savefig('%s/lnRatios.png'%directory, 
                              dpi=4*96)
    matplotlib.pyplot.close()
    
    y = numpy.log(
                    (Fobss**2)/((Fcalcs**2)*numpy.exp(-dB_2D*(q2Ds**2)))
                 )
    sig_y = 2*sigFobss/numpy.absolute(Fobss)
     
    x = numpy.linspace(0.0, 0.28, 1000) #0 
    l = lnA - dB_rod*x
    
    uniqueQrods = cStar*numpy.arange(16) #22
    print uniqueQrods
    
    y_avgs = []
    sig_y_avgs = []
    for uniqueQrod in uniqueQrods:
        y_bin     = [y[i]     for i in range(len(y)) if qRods[i] == uniqueQrod]
        sig_y_bin = [sig_y[i] for i in range(len(y)) if qRods[i] == uniqueQrod]
        
        y_avg = numpy.average(y_bin)
        y_avgs.append(y_avg)
        sig_y_bin = numpy.asarray(sig_y_bin)
        sig_y_avg = numpy.sqrt((sig_y_bin**2).sum())/len(y_bin)
        sig_y_avgs.append(sig_y_avg)
        print y_avg, sig_y_avg        
    
        
    f = matplotlib.pyplot.figure() #figsize=(20, 10))
    matplotlib.pyplot.scatter(qRods**2, y, s=2)
    matplotlib.pyplot.scatter(uniqueQrods**2, 
                              y_avgs,
                              facecolors='m', 
                              edgecolors='none', 
                              s=40, 
                              alpha=0.4)
    matplotlib.pyplot.errorbar(uniqueQrods**2, 
                               y_avgs,
                               yerr = sig_y_avgs,
                               linestyle="None",
                               color='m')
    matplotlib.pyplot.plot(x, l)
    matplotlib.pyplot.tick_params(axis='both', which='major', labelsize=20)
    myAxes = matplotlib.pyplot.gca()
    myAxes.set_xlim([-0.01,0.30]) # to comment
    myAxes.set_xlabel(r'q$^{2}$$_{\rm rod}$ [A$^{-2}$]', fontsize=22) 
    myAxes.set_ylabel(r'ln[I$_{\rm obs}$/(I$_{\rm calc}$exp(-$\Delta$B$_{\rm 2D}$q$^{2}$$_{\rm 2D}$))]', fontsize=22) 
    f.tight_layout()
    matplotlib.pyplot.savefig('%s/scatter.pdf'%directory, dpi=96*2)
    matplotlib.pyplot.savefig('%s/scatter.png'%directory, dpi=96*2)
    matplotlib.pyplot.close()
    

def anisotropy_noSigF(directory, inFile, cellSize, cAxis):
    
    cStar     = (2*numpy.pi)/cAxis
    
    qRods    = []
    q2Ds     = []
    Fobss    = []
    Fcalcs   = []
    ratios   = []
    
    inputFile = ('%s/%s'%(directory, inFile))
    inputFile_read = open(inputFile, 'r')
    for line in inputFile_read:
        splitLine = line.split()
        
        h       = int(splitLine[0])
        k       = int(splitLine[1])
        l       = int(splitLine[2])
        Fobs    = float(splitLine[3])
        Fcalc   = float(splitLine[4]) # PDB MODEL AND BULK SOLVENT, NO ANISOTROPY CORRECTION
        ratio   = (Fobs**2)/(Fcalc**2)
        
        if h == 0 and k == 0:
            continue
        
        qRod = l*cStar
        d2D  = simulate_resolution.resolution(cellSize, h, k, 0.0)
        q2D  = (2*numpy.pi)/d2D
        
        qRods.append(abs(qRod))
        q2Ds.append(q2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
        ratios.append(ratio)
           
    qRods    = numpy.asarray(qRods)
    q2Ds     = numpy.asarray(q2Ds)
    Fobss    = numpy.asarray(Fobss)
    Fcalcs   = numpy.asarray(Fcalcs)
    ratios   = numpy.asarray(ratios)

    lnRatios     = numpy.log(ratios)
    q2D_squared  = q2Ds**2
    qRod_squared = qRods**2
    
    popt_lnRatios, \
    pcov_lnRatios = scipy.optimize.curve_fit(lnRatioModel, 
                                             (q2D_squared, qRod_squared), 
                                             lnRatios)
                                            
    lnA    = popt_lnRatios[0]
    dB_2D  = popt_lnRatios[1]
    dB_rod = popt_lnRatios[2]
   
    print 'Parameters: ', lnA, dB_2D, dB_rod
    
    xx, yy = numpy.meshgrid(numpy.linspace(0, max(q2D_squared),  101),
                            numpy.linspace(0, max(qRod_squared), 101))
                            
    zz_ratio = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_ratio[i,j] = lnRatioModel((xx[i,j], yy[i,j]), 
                                         lnA, 
                                         dB_2D, 
                                         dB_rod)
    
    # 3D Figure
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(q2D_squared, 
                   qRod_squared, 
                   lnRatios, 
                   c='k', 
                   marker='o', 
                   s=0.5)
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
    myAxes.set_zlabel(r'ln(I$_{\rm obs}$/I$_{\rm calc}$)') 
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
    matplotlib.pyplot.tight_layout()  
    print 'Showing 3D figure.'
    matplotlib.pyplot.savefig('%s/lnRatios.png'%directory, 
                              dpi=4*96)
    matplotlib.pyplot.close()
    
    y = numpy.log(
                    (Fobss**2)/((Fcalcs**2)*numpy.exp(-dB_2D*(q2Ds**2)))
                 )
     
    x = numpy.linspace(0, 1.0, 1000) #0.28
    l = lnA - dB_rod*x
    
    uniqueQrods = cStar*numpy.arange(16) #16
    print uniqueQrods
    
    y_avgs = []
    for uniqueQrod in uniqueQrods:
        y_bin     = [y[i]     for i in range(len(y)) if qRods[i] == uniqueQrod]
        
        y_avg = numpy.average(y_bin)
        y_avgs.append(y_avg)
        
        print y_avg     
    
        
    matplotlib.pyplot.figure() #figsize=(20, 10))
    matplotlib.pyplot.scatter(qRods**2, y, s=2)
    matplotlib.pyplot.scatter(uniqueQrods**2, 
                              y_avgs,
                              facecolors='m', 
                              edgecolors='none', 
                              s=40, 
                              alpha=0.4)
    matplotlib.pyplot.plot(x, l)
    matplotlib.pyplot.tick_params(axis='both', which='major', labelsize=15)
    myAxes = matplotlib.pyplot.gca()
    myAxes.set_xlim([-0.01,0.30]) # TO COMMENT....
    myAxes.set_xlabel(r'q$^{2}$$_{\rm rod}$ [A$^{-2}$]', fontsize=25) 
    myAxes.set_ylabel(r'ln[I$_{\rm obs}$/(I$_{\rm calc}$exp(-$\Delta$B$_{\rm 2D}$q$^{2}$$_{\rm 2D}$))]', fontsize=25) 
    matplotlib.pyplot.savefig('%s/Scatter.png'%directory, dpi=96*2)
    matplotlib.pyplot.close()  
    
 
    
if __name__ == "__main__":
    
    print "\n**** CALLING anisotropyCheck ****"  

###2D-SFX###    
    directory  = './Anisotropy'
    inputFile  = 'model_no_anisotropicCorr_refine_001_f_model_ed.txt'
    outputFile = 'f_model_onlyBulkSolvent_B_from_5b6v.txt'
    cellSize   = 62.45     # A
    cAxis      = 180       # Double sampling!

###5B6V###    
#    directory  = './Anisotropy_5b6v'
#    inputFile  = '5b6v_edited_refine_001_f_model_edited.txt'
#    outputFile = 'f_model_onlyBulkSolvent_5b6v.txt'
#    cellSize   = 62.50     # A
#    cAxis      = 112       # A
    
    #prepareFile(directory, inputFile, outputFile)    
    anisotropy(directory, outputFile, cellSize, cAxis)

###1BRD###    
#    directory  = './Anisotropy_1brd'
#    inputFile  = '1brd_refine_001_f_model_edited.txt'
#    outputFile = 'f_model_onlyBulkSolvent_1brd.txt'
#    cellSize   = 62.45     # A
#    cAxis      = 100       # A
    
    #prepareFile_noSigF(directory, inputFile, outputFile)
    #anisotropy_noSigF(directory, outputFile, cellSize, cAxis)