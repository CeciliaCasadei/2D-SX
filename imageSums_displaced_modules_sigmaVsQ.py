# -*- coding: utf-8 -*-
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import scipy.optimize
import pickle
import sys
import getopt

import imageSums_utilities



def sigmaVsQ(myArguments):

    str1 = "--selectedRun <selectedRun> --photonThreshold <photonThreshold> --dThreshold <dThreshold>"
    str2 = "--nCountsPerPhoton <nCountsPerPhoton> --cellSize <cellSize>"
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "photonThreshold=",
                                                                 "dThreshold=",
                                                                 "nCountsPerPhoton=",
                                                                 "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules_sigmaVsQ.py %s %s'%(str1, str2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules_sigmaVsQ.py %s %s'%(str1, str2)
            sys.exit()
        elif option == "--selectedRun":
            runNumber = value.zfill(4)
        elif option == "--photonThreshold":
            photonThreshold = float(value)
        elif option == "--dThreshold":
            dThreshold = float(value)
        elif option == "--nCountsPerPhoton":
            nCountsPerPhoton = int(value)
        elif option == "--cellSize":
            cellSize = float(value)
            
        
    
    # TO BE FITTED
    qs = []
    sigXs= []
    sigYs= []
    
    # FOLDERS
    outputFolder = './Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits'%runNumber
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
    
    # RECIPROCAL CELL
    
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I                                                  # A^(-1)
    
    # LOOP ON SPOTS WITH GAUSS INTENSITY ABOVE 1 PHOTON AND GOOD GAUSSIAN FIT
    fOpen = open('./Output_r%s/Output_imageSums_moduleDisplacements/imageSums_displaced_modules.txt'%runNumber, 'r')
    for hk_line in fOpen:
        splittedLine = hk_line.split()
        try:
            h = int(splittedLine[0])
        except:
            continue
        
        h              = int(splittedLine[0])
        k              = int(splittedLine[1])
        gaussAmplitude = float(splittedLine[2])
        Dx0            = float(splittedLine[3])
        Dy0            = float(splittedLine[4])    
        sigX           = float(splittedLine[5])
        sigY           = float(splittedLine[6])
        Igauss         = float(splittedLine[7])
        Isum           = float(splittedLine[8])
        
        reciprocalVector = [h, k]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
        
        if gaussAmplitude < 0:
            print '(%4d, %4d) DISCARDED: negative Gauss amplitude'%(h, k)
            continue    
        if abs(Igauss) < photonThreshold*nCountsPerPhoton:
            print '(%4d, %4d) DISCARDED: I < %f photons'%(h, k, photonThreshold)
            continue
        if Isum < photonThreshold:
            print '(%4d, %4d) DISCARDED: I < %f photons'%(h, k, photonThreshold)
            continue
        if sigX > 5 or sigY > 5:
            print '(%4d, %4d) DISCARDED: Gauss width is too large.'%(h, k)
            continue
        if q < 0.15:
            continue
        if sigX < 0:
            print 'NEGATIVE sigma_x', h, k, q, sigX
            continue
        if sigY < 0:
            print 'NEGATIVE sigma_y', h, k, q, sigY
            continue
        if (Dx0 > dThreshold or Dy0 > dThreshold):
            print 'Large distance:', h, k, q, Dx0, Dy0
            continue
        
        qs.append(q)
        sigXs.append(sigX)
        sigYs.append(sigY)
        print 'KEEPING: h=%d k=%d q=%.2f sigX=%.1f sigY=%.1f'%(h, k, q, sigX, sigY)
    
       
    ### PLOT SIGMA VS Q ###
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.savefig('%s/sigXvsQ_%.1fph.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.savefig('%s/sigYvsQ_%.1fph.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    
    
    # FIT sigX VS q WITH A 4th ORDER POLY
    popt_sigX, pcov = scipy.optimize.curve_fit(imageSums_utilities.poly_4, qs, sigXs)
    x = numpy.linspace(0.1, 1.6, 100)
    y_poly_4 = imageSums_utilities.poly_4(x, *popt_sigX)
    
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.title('%s'%popt_sigX)
    matplotlib.pyplot.plot(x, y_poly_4, label='fit')
    matplotlib.pyplot.savefig('%s/sigXvsQ_poly_4_%.1fph.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    
    
    # FIT sigY VS q WITH A 4th ORDER POLY
    popt_sigY, pcov = scipy.optimize.curve_fit(imageSums_utilities.poly_4, qs, sigYs)
    x = numpy.linspace(0.1, 1.6, 100)
    y_poly_4 = imageSums_utilities.poly_4(x, *popt_sigY)
    
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.title('%s'%popt_sigY)
    matplotlib.pyplot.plot(x, y_poly_4, label='fit')
    matplotlib.pyplot.savefig('%s/sigYvsQ_poly_4_%.1fph.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    offset_sigX = popt_sigX[0]
    offset_sigY = popt_sigY[0]
    
    avg_offset = (offset_sigX + offset_sigY)/2
    print 'AVERAGE OFFSET: ', avg_offset
    
    # FIT sigX VS q WITH A 4th ORDER POLY, FIXED OFFSET
    popt_sigX, pcov = scipy.optimize.curve_fit(lambda (x), a, b: imageSums_utilities.poly_4((x), avg_offset, a, b), (qs), sigXs) 
    x = numpy.linspace(0.1, 1.6, 100)
    y_poly_4_sigX = imageSums_utilities.poly_4(x, avg_offset, popt_sigX[0], popt_sigX[1])
    
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigXs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.title('%s'%popt_sigX)
    matplotlib.pyplot.plot(x, y_poly_4_sigX, label='fit')
    matplotlib.pyplot.savefig('%s/sigXvsQ_poly_4_%.1fph_fixed_offset.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    # FIT sigY VS q WITH A 4th ORDER POLY, FIXED OFFSET
    popt_sigY, pcov = scipy.optimize.curve_fit(lambda (x), a, b: imageSums_utilities.poly_4((x), avg_offset, a, b), (qs), sigYs) 
    x = numpy.linspace(0.1, 1.6, 100)
    y_poly_4_sigY = imageSums_utilities.poly_4(x, avg_offset, popt_sigY[0], popt_sigY[1])
    
    matplotlib.pyplot.figure(figsize = (20, 5))
    matplotlib.pyplot.scatter(qs, sigYs, color='b', label='data', alpha=0.5)
    axes = matplotlib.pyplot.gca()
    axes.set_xlim([0, 1.8])
    axes.set_ylim([0.5, 4.0])
    matplotlib.pyplot.title('%s'%popt_sigY)
    matplotlib.pyplot.plot(x, y_poly_4_sigY, label='fit')
    matplotlib.pyplot.savefig('%s/sigYvsQ_poly_4_%.1fph_fixed_offset.png'%(outputFolder, photonThreshold))
    matplotlib.pyplot.close()
    
    ### SUMMARY PLOT ### 
    myFigure, (ax1, ax2) = matplotlib.pyplot.subplots(2, 1)    
    ax1.scatter(qs, sigXs  , color='b', label='data', alpha=0.5, s=3) 
    ax1.set_ylabel(r'$\sigma_{\rm rad}$ (pix)', fontsize = 17) 
    ax1.set_ylim([0.5, 4.0]) 
    ax1.text(0.02, 0.96, '(a)', 
             horizontalalignment='left', verticalalignment='top', fontsize=25, transform = ax1.transAxes) 
    ax1.text(0.96, 0.95, r'$\sigma_{\rm rad} = %.2f + %.2f q^2 + %.2f q^4$'%(avg_offset, popt_sigX[0], popt_sigX[1]), 
             horizontalalignment='right', verticalalignment='top', fontsize=15, transform = ax1.transAxes) 
    ax1.tick_params(axis='both', which='major', labelsize=12, length=6)
    ax1.plot(x, y_poly_4_sigX, label='fit', color='m') 
    
    ax2.scatter(qs, sigYs, color='b', label='data', alpha=0.5, s=3) 
    ax2.set_xlabel(r'$q$ ($\AA^{-1}$)', fontsize = 17) 
    ax2.set_ylabel(r'$\sigma_{\rm azi}$ (pix)', fontsize = 17) 
    ax2.set_ylim([0.5, 4.0]) 
    ax2.text(0.02, 0.96, '(b)', 
             horizontalalignment='left', verticalalignment='top', fontsize=25, transform = ax2.transAxes) 
    ax2.text(0.96, 0.95, r'$\sigma_{\rm azi} = %.2f + %.2f q^2 + %.2f q^4$'%(avg_offset, popt_sigY[0], popt_sigY[1]), 
             horizontalalignment='right', verticalalignment='top', fontsize=15, transform = ax2.transAxes) 
    ax2.tick_params(axis='both', which='major', labelsize=12, length=6) 
    ax2.plot(x, y_poly_4_sigY, label='fit', color='m') 
    
    matplotlib.pyplot.tight_layout()                        
    matplotlib.pyplot.savefig('%s/fits_poly_4_fixed_offset_threshold_%.1f_ph_dThresh_%.1f.png'%(outputFolder, photonThreshold, dThreshold), dpi=4*96) 
    matplotlib.pyplot.savefig('%s/fits_poly_4_fixed_offset_threshold_%.1f_ph_dThresh_%.1f.pdf'%(outputFolder, photonThreshold, dThreshold), dpi=4*96) 
    matplotlib.pyplot.close()
    
    # SAVE SIGMA FIT PARAMETERS
    sigmaX_curve = open('%s/sigmaXCurveParameters.pkl'%outputFolder, 'wb')
    sigX_parameters = [avg_offset, popt_sigX[0], popt_sigX[1]]
    pickle.dump(sigX_parameters, sigmaX_curve)
    sigmaX_curve.close()
    
    sigmaY_curve = open('%s/sigmaYCurveParameters.pkl'%outputFolder, 'wb')
    sigY_parameters = [avg_offset, popt_sigY[0], popt_sigY[1]]
    pickle.dump(sigY_parameters, sigmaY_curve)
    sigmaY_curve.close()

if __name__ == "__main__":
    print "\n**** CALLING imageSums_displaced_modules_sigmaVsQ ****"
    sigmaVsQ(sys.argv[1:])   