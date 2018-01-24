# -*- coding: utf-8 -*-
import numpy
import pickle
import sys
import getopt
import scipy.optimize
import matplotlib.pyplot

import imageSums_utilities


            
def sumTilted_Function_spotWidthModeling(myArguments):
    nCountsPerPhoton = 26
    
    # CRITERIA
    dMax = 3
    sigmaMax = 4
    sigmaMin = 0.8
    gaussImin = 1
    nTerms_min = 20
    
    str1 = '--halfWidth <halfWidth> --tiltAngle <tiltAngle> --nominalCell <nominalCell>'
           
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["halfWidth=", 
                                                                 "tiltAngle=",
                                                                 "nominalCell="
                                                                 ])
    except getopt.GetoptError:
        print ('Error Usage: python imageSumming_tilted_spotWidthModeling.py %s'
               %(str1))
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print ('Usage: python imageSumming_tilted_spotWidthModeling.py %s'
                   %(str1))
            sys.exit()
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--tiltAngle":
            tiltAngle_deg = value
        elif option == "--nominalCell":
            nominalCell = float(value)
  
    # FOLDERS  
    inputFolder = './Output_imageSum_gaussFit_tilt_%s'%tiltAngle_deg  

    # RECIPROCAL CELL
    directCell = nominalCell * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],
                                             [0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I  
    
    # LOAD PEAK SUMS    
    sums_file = open('%s/sumMatrix_dictionary_list_gaussFit_tilt_%s.pkl'
                     %(inputFolder, tiltAngle_deg), 'rb')
    dictionaryList = pickle.load(sums_file)
    sums_file.close()    
    print '%d peak sums.'%len(dictionaryList)
    
    # EXTRACT sigma_x, sigma_y
    sigmaXs = []
    sigmaYs = []
    qRods = []
    q2Ds = []    
    for dictionary in dictionaryList:
        h = dictionary['h']
        k = dictionary['k']
        qRod = dictionary['qRod']
        nTerms = dictionary['nTerms']
        sigmaX = dictionary['refined_sigma_x'] 
        sigmaY = dictionary['refined_sigma_y'] 
        x0 = dictionary['refined_x0'] 
        y0 = dictionary['refined_y0'] 
        amplitude = dictionary['refined_amplitude'] 
        gaussIntegral = dictionary['gauss_integral']
        
        if numpy.isnan(gaussIntegral):
            continue
        
        d = numpy.sqrt((x0 - halfWidth)**2 + (y0 - halfWidth)**2)
        
        # CRITERIA
        if (sigmaX < sigmaMax and
            sigmaY < sigmaMax and
            sigmaX > sigmaMin and
            sigmaY > sigmaMin and
            d < dMax and
            gaussIntegral > nCountsPerPhoton*gaussImin and
            amplitude > 0 and
            nTerms > nTerms_min):
                          
            reciprocalVector = [h, k]*reciprocalCellRows
            q_x = reciprocalVector[0,0]         # A^(-1)
            q_y = reciprocalVector[0,1]         # A^(-1)
            q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1) 
            
            sigmaXs.append(sigmaX)
            sigmaYs.append(sigmaY)
            qRods.append(qRod)
            q2Ds.append(q_2D)
        
    sigmaXs = numpy.asarray(sigmaXs)
    sigmaYs = numpy.asarray(sigmaYs)
    qRods   = numpy.asarray(qRods)
    q2Ds    = numpy.asarray(q2Ds)
    
    # PLOT sigma_x, sigma_y as a function of 2_2D, q_z.
    matplotlib.pyplot.figure()  
    matplotlib.pyplot.title(r'$\sigma_x$', fontsize = 22)      
    matplotlib.pyplot.scatter(q2Ds, qRods, c=sigmaXs, 
                              vmin=1, vmax=4, 
                              s=25, edgecolors='none')
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.2])
    axes.set_ylim([-0.5, +0.5])
    matplotlib.pyplot.cool()
    matplotlib.pyplot.colorbar()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.savefig('%s/sigmaX_tilt_%s'%(inputFolder, tiltAngle_deg))
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.title(r'$\sigma_y$', fontsize = 22)          
    matplotlib.pyplot.scatter(q2Ds, qRods, c=sigmaYs, 
                              vmin=1, vmax=4, 
                              s=25, edgecolors='none')
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.2])
    axes.set_ylim([-0.5, +0.5])
    matplotlib.pyplot.cool()
    matplotlib.pyplot.colorbar()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.savefig('%s/sigmaY_tilt_%s'%(inputFolder, tiltAngle_deg))
    matplotlib.pyplot.close()
    
    # MODEL sigma_x = offset + a*q_2D**2 + b*q_z**2
    initial_guess = (0.8, 0, 0)
    popt_sigX, \
    pcov_sigX = scipy.optimize.curve_fit(imageSums_utilities.poly_2_2D, 
                                         (q2Ds, qRods), 
                                          sigmaXs, 
                                          p0=initial_guess)
                                          
    xx, yy = numpy.meshgrid(numpy.linspace(0,1.2, 1001),
                            numpy.linspace(-0.5,0.5, 1001))

    zz_sigX = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_sigX[i,j] = imageSums_utilities.poly_2_2D((xx[i,j], yy[i,j]), 
                                                          popt_sigX[0], 
                                                          popt_sigX[1], 
                                                          popt_sigX[2])  
    
    # PLOT THE MODELING FUNCTION
    matplotlib.pyplot.title(r'$\sigma_x$ = %.2f + %.2fq$_{\rm 2D}$$^2$ +%.2fq$_{z}$$^2$'
                            %(popt_sigX[0], popt_sigX[1], popt_sigX[2]))
    matplotlib.pyplot.pcolor(xx, yy, zz_sigX, vmin=1, vmax=4)
    matplotlib.pyplot.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    matplotlib.pyplot.colorbar()
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=5, edgecolors='none')    
    matplotlib.pyplot.savefig('%s/fit_sigmaX_tilt_%s'%(inputFolder, tiltAngle_deg), 
                              dpi=96*4)
    matplotlib.pyplot.close()
    
    # MODEL sigma_y = offset + a*q_2D**2 + b*q_z**2
    popt_sigY, \
    pcov_sigY = scipy.optimize.curve_fit(imageSums_utilities.poly_2_2D, 
                                         (q2Ds, qRods), 
                                          sigmaYs, 
                                          p0=initial_guess)
   
    zz_sigY = numpy.zeros(xx.shape)
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            zz_sigY[i,j] = imageSums_utilities.poly_2_2D((xx[i,j], yy[i,j]), 
                                                          popt_sigY[0], 
                                                          popt_sigY[1], 
                                                          popt_sigY[2])  
    
    # PLOT THE MODELING FUNCTION
    matplotlib.pyplot.title(r'$\sigma_y$ = %.2f + %.2fq$_{\rm 2D}$$^2$ +%.2fq$_{z}$$^2$'
                            %(popt_sigY[0], popt_sigY[1], popt_sigY[2]))
    matplotlib.pyplot.pcolor(xx, yy, zz_sigY, vmin=1, vmax=4)
    matplotlib.pyplot.axis([xx.min(), xx.max(), yy.min(), yy.max()])
    matplotlib.pyplot.colorbar()
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
    axes.set_ylabel(r'q$_z$ [A$^{-1}$]', fontsize = 17) 
    matplotlib.pyplot.scatter(q2Ds, qRods, c='k', s=5, edgecolors='none')    
    matplotlib.pyplot.savefig('%s/fit_sigmaY_tilt_%s'%(inputFolder, tiltAngle_deg), 
                              dpi=96*4)
    matplotlib.pyplot.close()   
    
    # SAVE SIGMA FIT PARAMETERS
    sigmaX_curve = open('%s/sigmaXCurveParameters.pkl'%inputFolder, 'wb')
    sigX_parameters = [popt_sigX[0], popt_sigX[1], popt_sigX[2]]
    pickle.dump(sigX_parameters, sigmaX_curve)
    sigmaX_curve.close()
    
    sigmaY_curve = open('%s/sigmaYCurveParameters.pkl'%inputFolder, 'wb')
    sigY_parameters = [popt_sigY[0], popt_sigY[1], popt_sigY[2]]
    pickle.dump(sigY_parameters, sigmaY_curve)
    sigmaY_curve.close()
    
    # PLOT SLICES AT FIXED qRod 
    qRodBins = numpy.arange(-0.5, +0.5, 0.01)
    for qRodBin in qRodBins:
        qRodValues = [qRods[i] for i in range(0, len(qRods)) 
                               if abs(qRods[i] - qRodBin) <= 0.001]
        if len(qRodValues) > 0:
            q2DValues  = [q2Ds[i]    for i in range(0, len(qRods)) 
                                     if abs(qRods[i] - qRodBin) <= 0.001]
            sigXvalues = [sigmaXs[i] for i in range(0, len(qRods)) 
                                     if abs(qRods[i] - qRodBin) <= 0.001]
            sigYvalues = [sigmaYs[i] for i in range(0, len(qRods)) 
                                     if abs(qRods[i] - qRodBin) <= 0.001]
            print 'qRod = %.2f, %d data points'%(qRodBin, len(sigXvalues))
            
            indices_list = numpy.argwhere(abs(yy-qRodBin) <= 0.00075)
            print 'Modeling function calculated on %d points'%len(indices_list)
            
            q2DValues_fit  = []
            sigXvalues_fit = []
            sigYvalues_fit = []
            for indices in indices_list:
                q2DValues_fit.append(xx[indices[0], indices[1]])
                sigXvalues_fit.append(zz_sigX[indices[0], indices[1]]) 
                sigYvalues_fit.append(zz_sigY[indices[0], indices[1]])   
                               
            matplotlib.pyplot.figure()
            matplotlib.pyplot.title(r'q$_z$=%.2f A$^{-1}$'%qRodBin)
            matplotlib.pyplot.scatter(q2DValues, sigXvalues, 
                                      s=20, edgecolors='none', alpha=0.5)
            matplotlib.pyplot.scatter(q2DValues_fit, sigXvalues_fit, 
                                      s=5, edgecolors='none', c='c')
            axes = matplotlib.pyplot.gca()
            axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
            axes.set_ylabel(r'$\sigma_x$ [pxls]', fontsize = 17) 
            axes.set_xlim([0,   1.2])
            axes.set_ylim([1.0, 4.0])
            matplotlib.pyplot.savefig('%s/fit_sigmaX_tilt_%s_qRod_%.2f.png'
                                       %(inputFolder, tiltAngle_deg, qRodBin), 
                                       dpi=96*4)
            matplotlib.pyplot.close()
            
            matplotlib.pyplot.figure()
            matplotlib.pyplot.title(r'q$_z$=%.2f A$^{-1}$'%qRodBin)
            matplotlib.pyplot.scatter(q2DValues, sigYvalues, 
                                      s=20, edgecolors='none', alpha=0.5)
            matplotlib.pyplot.scatter(q2DValues_fit, sigYvalues_fit, 
                                      s=5, edgecolors='none', c='c')          
            axes = matplotlib.pyplot.gca()
            axes.set_xlabel(r'q$_{\rm 2D}$ [A$^{-1}$]', fontsize = 17) 
            axes.set_ylabel(r'$\sigma_y$ [pxls]', fontsize = 17) 
            axes.set_xlim([0,   1.2])
            axes.set_ylim([1.0, 4.0])
            matplotlib.pyplot.savefig('%s/fit_sigmaY_tilt_%s_qRod_%.2f.png'
                                       %(inputFolder, tiltAngle_deg, qRodBin), 
                                       dpi=96*4)
            matplotlib.pyplot.close()                        
   
        

if __name__ == "__main__":
    print "\n**** CALLING imageSumming_tilted_spotWidthModeling ****"
    sumTilted_Function_spotWidthModeling(sys.argv[1:])   