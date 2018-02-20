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
        
    print len(qRods), len(q2Ds), len(Fobss), len(Fcalcs)
    print max(qRods), max(q2Ds)
    print max(Fobss), max(Fcalcs)
    print numpy.average(numpy.asarray(Fobss)), numpy.average(numpy.asarray(Fcalcs))
    
    
    grid_x, \
    grid_y = numpy.meshgrid( numpy.linspace(-0.99*max(q2Ds),  0.99*max(q2Ds),  100), 
                             numpy.linspace(-0.99*max(qRods), 0.99*max(qRods), 100))

    print grid_x
    print grid_y
    
    xy = [q2Ds,qRods]
    xy = numpy.asarray(xy)
    xy = xy.T
    
    print xy.shape

    grid_Fobss  = scipy.interpolate.griddata(xy, Fobss,  (grid_x, grid_y), method='linear')
    grid_Fcalcs = scipy.interpolate.griddata(xy, Fcalcs, (grid_x, grid_y), method='linear')
    
    print grid_x.shape, grid_y.shape, grid_Fobss.shape, grid_Fcalcs.shape

    guess = (max(Fobss), 1.0, 1.0)    
    popt, \
    pcov = scipy.optimize.curve_fit(twoD_Gaussian_simple, 
                                   (grid_x, grid_y), 
                                    grid_Fobss.ravel(), 
                                    p0=guess)
    print popt
    A_obs = popt[0]
    sig_2D_obs = popt[1]
    sig_rod_obs = popt[2]
    print sig_2D_obs/sig_rod_obs
    
    
    guess = (max(Fcalcs), 1.0, 1.0)
    
    popt, \
    pcov = scipy.optimize.curve_fit(twoD_Gaussian_simple, 
                                   (grid_x, grid_y), 
                                    grid_Fcalcs.ravel(), 
                                    p0=guess)
    print popt    
    A_calc = popt[0]
    sig_2D_calc = popt[1]
    sig_rod_calc = popt[2]
    print sig_2D_calc/sig_rod_calc
    
    print min(q2Ds), max(q2Ds)
    print min(qRods), max(qRods)
    
    q2D_bins = numpy.linspace(0, max(q2Ds), 50)
    for q2D_bin_n in range(0, len(q2D_bins)-1):
        
        left  = q2D_bins[q2D_bin_n]
        right = q2D_bins[q2D_bin_n+1]
        middle = 0.5*(left+right)
        
        qrod_bin  = [abs(qRods[i])  for i in range(0, len(q2Ds)) 
                                    if left < abs(q2Ds[i]) < right]
        Fobs_bin  = [Fobss[i]       for i in range(0, len(q2Ds)) 
                                    if left < abs(q2Ds[i]) < right]
        Fcalc_bin = [Fcalcs[i]      for i in range(0, len(q2Ds)) 
                                    if left < abs(q2Ds[i]) < right]
        
        if len(Fobs_bin) > 0:
            
            matplotlib.pyplot.figure()
            matplotlib.pyplot.scatter(qrod_bin, Fobs_bin,  color='b')
            matplotlib.pyplot.scatter(qrod_bin, Fcalc_bin, color='r')
            
            # Plot fit:
            y = numpy.linspace(min(qrod_bin), max(qrod_bin), 1000)
            x0 = middle*numpy.ones(y.shape)
            obs  = twoD_Gaussian_simple((x0, y), A_obs,  sig_2D_obs,  sig_rod_obs)
            calc = twoD_Gaussian_simple((x0, y), A_calc, sig_2D_calc, sig_rod_calc)
            matplotlib.pyplot.plot(y, obs,  'b')
            matplotlib.pyplot.plot(y, calc, 'r')           
            
            axes = matplotlib.pyplot.gca()
            axes.set_xlabel(r'q$_{\rm rod}$ [A$^{-1}$]', fontsize = 17) 
            axes.set_ylabel(r'F', fontsize = 17) 
            matplotlib.pyplot.title('q2D bin: %.2f - %.2f A-1'%(left, right))
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig('%s/FvsQrod_fixed_q2D_%d'%(directory,
                                                                 q2D_bin_n))
            matplotlib.pyplot.close()
           
    qRod_bins = numpy.linspace(0, max(qRods), 70)
    for qRod_bin_n in range(0, len(qRod_bins)-1):
        
        left  = qRod_bins[qRod_bin_n]
        right = qRod_bins[qRod_bin_n+1]
        middle = 0.5*(left+right)
        
        q2Ds_bin  = [abs(q2Ds[i]) for i in range(0, len(q2Ds)) 
                                  if left < abs(qRods[i]) < right]
        Fobs_bin  = [Fobss[i]     for i in range(0, len(q2Ds)) 
                                  if left < abs(qRods[i]) < right]
        Fcalc_bin = [Fcalcs[i]    for i in range(0, len(q2Ds)) 
                                  if left < abs(qRods[i]) < right]
        
        if len(Fobs_bin) > 0:
            
            matplotlib.pyplot.figure()
            matplotlib.pyplot.scatter(q2Ds_bin, Fobs_bin,  color='b')
            matplotlib.pyplot.scatter(q2Ds_bin, Fcalc_bin, color='r')
            
            # Plot fit:
            x = numpy.linspace(min(q2Ds_bin), max(q2Ds_bin), 1000)
            y0 = middle*numpy.ones(x.shape)
            obs  = twoD_Gaussian_simple((x, y0), A_obs,  sig_2D_obs,  sig_rod_obs)
            calc = twoD_Gaussian_simple((x, y0), A_calc, sig_2D_calc, sig_rod_calc)
            matplotlib.pyplot.plot(x, obs,  'b')
            matplotlib.pyplot.plot(x, calc, 'r')
            
            axes = matplotlib.pyplot.gca()
            axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
            axes.set_ylabel(r'F', fontsize = 17) 
            matplotlib.pyplot.title('qRod bin: %.2f - %.2f A-1'%(left, right))
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig('%s/Fvsq2D_fixed_qRod_%d'%(directory,
                                                                 qRod_bin_n))
            matplotlib.pyplot.close()


    
if __name__ == "__main__":
    print "\n**** CALLING anisotropyCheck ****"    
    anisotropy()
    
