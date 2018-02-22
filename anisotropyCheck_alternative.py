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
    
def ratioModel((x, y), amplitude, B_2D, B_z):  
    g = amplitude * numpy.exp( - (B_2D*x**2)
                               - (B_z *y**2)  )
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
    
    inputFile = ('%s/f_model.txt'%directory)
    inputFile_read = open(inputFile, 'r')
    for line in inputFile_read:
        splitLine = line.split()
        h     = int(splitLine[0])
        k     = int(splitLine[1])
        l     = int(splitLine[2])
        Fobs  = float(splitLine[3])
        Fcalc = float(splitLine[4])
        ratio = (Fobs**2)/(Fcalc**2)
        #ratio = Fobs/Fcalc
        
        qRod = l*cStar
        d_2D = simulate_resolution.resolution(cellSize, h, k, 0.0)
        q_2D = (2*numpy.pi)/d_2D
        
        qRods.append(qRod)
        q2Ds.append(q_2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
        ratios.append(ratio)
     
    # For easier gaussian fitting.
    for i in range(0, len(Fobss)):
        qRod = qRods[i]
        q2D = - q2Ds[i]
        Fobs = Fobss[i]
        Fcalc = Fcalcs[i]
        ratio = ratios[i]
        
        qRods.append(qRod)
        q2Ds.append(q2D)
        Fobss.append(Fobs)
        Fcalcs.append(Fcalc)
        ratios.append(ratio)
        
    Fobss   = numpy.asarray(Fobss)
    Fcalcs  = numpy.asarray(Fcalcs)
    qRods   = numpy.asarray(qRods)
    q2Ds    = numpy.asarray(q2Ds)
    ratios  = numpy.asarray(ratios)
    
    # PLOT sigma_x, sigma_y as a function of q_2D, q_z.
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
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=1, edgecolors='none')
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
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=1, edgecolors='none') 
    axes.set_aspect('equal')
    fig.canvas.draw()
    matplotlib.pyplot.tight_layout()    
    matplotlib.pyplot.savefig('%s/Fit_Fcalc.png'%(directory), 
                              dpi=96*4)
    matplotlib.pyplot.close()   


    #MODEL 
    initial_guess = (max(ratios), 1.0, 1.0)
    popt_ratios, \
    pcov_ratios = scipy.optimize.curve_fit(ratioModel, 
                                          (q2Ds, qRods), 
                                           ratios, 
                                           p0=initial_guess)
   
    zz_ratio = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_ratio[i,j] = ratioModel((xx[i,j], yy[i,j]), 
                                        popt_ratios[0], 
                                        popt_ratios[1], 
                                        popt_ratios[2])  
                                                          
    print popt_ratios
    
    # PLOT THE MODELING FUNCTION
    fig, axes = matplotlib.pyplot.subplots(figsize=(5,5), facecolor='c')
    matplotlib.pyplot.title(r'Ratios')
    matplotlib.pyplot.pcolor(xx, yy, zz_ratio)
    matplotlib.pyplot.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    matplotlib.pyplot.colorbar()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=1, edgecolors='none') 
    axes.set_aspect('equal')
    fig.canvas.draw()
    matplotlib.pyplot.tight_layout()    
    matplotlib.pyplot.savefig('%s/Ratios.png'%(directory), 
                              dpi=96*4)
    matplotlib.pyplot.close()   
    
    
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(q2Ds, qRods, ratios, c='r', marker='o', s=6)
    mySurface = myAxes.plot_surface(xx, yy, zz_ratio, rstride=1, cstride=1, cmap=cm.coolwarm,
                                    linewidth=0, antialiased=False, alpha = 0.3)
    myFigure.colorbar(mySurface, shrink=0.5, aspect=5)
  
    
    myAxes.zaxis.set_major_locator(LinearLocator(10))
    myAxes.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    myAxes.set_zlim(0, 40)
    myAxes.set_xlabel('X')
    myAxes.set_ylabel('Y')
    myAxes.set_zlabel('Z')  
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
   
    
    print 'show'
    matplotlib.pyplot.show()
    
#    q2D_bins = numpy.linspace(0, max(q2Ds), 50)
#    for q2D_bin_n in range(0, len(q2D_bins)-1):
#        
#        left  = q2D_bins[q2D_bin_n]
#        right = q2D_bins[q2D_bin_n+1]
#        middle = 0.5*(left+right)
#        
#        qrod_bin   = [abs(qRods[i])  for i in range(0, len(q2Ds)) 
#                                     if left < abs(q2Ds[i]) < right]
#        ratio_bin  = [ratios[i]      for i in range(0, len(q2Ds)) 
#                                     if left < abs(q2Ds[i]) < right]
#        
#        if len(ratio_bin) > 0:
#            
#            matplotlib.pyplot.figure()
#            matplotlib.pyplot.scatter(qrod_bin, ratio_bin,  color='b')
#            
#            # Plot fit:
#            y = numpy.linspace(min(qrod_bin), max(qrod_bin), 1000)
#            x0 = middle*numpy.ones(y.shape)
#            ratio_fit  = ratioModel((x0, y), popt_ratios[0],  
#                                             popt_ratios[1],  
#                                             popt_ratios[2])
#            matplotlib.pyplot.plot(y, ratio_fit,  'b')
#                      
#            axes = matplotlib.pyplot.gca()
#            axes.set_xlabel(r'q$_{\rm rod}$ [A$^{-1}$]', fontsize = 17) 
#            axes.set_ylabel(r'Ratio', fontsize = 17) 
#            matplotlib.pyplot.title('q2D bin: %.2f - %.2f A-1'%(left, right))
#            matplotlib.pyplot.tight_layout()
#            matplotlib.pyplot.savefig('%s/ratioVsQrod_fixed_q2D_%d'%(directory,
#                                                                     q2D_bin_n))
#            matplotlib.pyplot.close()
           
#    qRod_bins = numpy.linspace(0, max(qRods), 70)
#    for qRod_bin_n in range(0, len(qRod_bins)-1):
#        
#        left  = qRod_bins[qRod_bin_n]
#        right = qRod_bins[qRod_bin_n+1]
#        middle = 0.5*(left+right)
#        
#        q2Ds_bin  = [abs(q2Ds[i]) for i in range(0, len(q2Ds)) 
#                                  if left < abs(qRods[i]) < right]
#        Fobs_bin  = [Fobss[i]     for i in range(0, len(q2Ds)) 
#                                  if left < abs(qRods[i]) < right]
#        Fcalc_bin = [Fcalcs[i]    for i in range(0, len(q2Ds)) 
#                                  if left < abs(qRods[i]) < right]
#        
#        if len(Fobs_bin) > 0:
#            
#            matplotlib.pyplot.figure()
#            matplotlib.pyplot.scatter(q2Ds_bin, Fobs_bin,  color='b')
#            matplotlib.pyplot.scatter(q2Ds_bin, Fcalc_bin, color='r')
#            
#            # Plot fit:
#            x = numpy.linspace(min(q2Ds_bin), max(q2Ds_bin), 1000)
#            y0 = middle*numpy.ones(x.shape)
#            ratio_fit  = twoD_Gaussian_simple((x, y0), A_obs,  sig_2D_obs,  sig_rod_obs)
#            calc = twoD_Gaussian_simple((x, y0), A_calc, sig_2D_calc, sig_rod_calc)
#            matplotlib.pyplot.plot(x, obs,  'b')
#            matplotlib.pyplot.plot(x, calc, 'r')
#            
#            axes = matplotlib.pyplot.gca()
#            axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
#            axes.set_ylabel(r'F', fontsize = 17) 
#            matplotlib.pyplot.title('qRod bin: %.2f - %.2f A-1'%(left, right))
#            matplotlib.pyplot.tight_layout()
#            matplotlib.pyplot.savefig('%s/Fvsq2D_fixed_qRod_%d'%(directory,
#                                                                 qRod_bin_n))
#            matplotlib.pyplot.close()

    
if __name__ == "__main__":
    print "\n**** CALLING anisotropyCheck ****"  
    #prepareFile()
    anisotropy()
    
