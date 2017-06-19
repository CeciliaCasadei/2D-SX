# -*- coding: utf-8 -*-
import math
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
from matplotlib import rcParams

### PAPER FIGURE (IMG N 74) ###
### NB Final orientation will be p -> LABEL INDICES NEED TO BE PERMUTED ###

def zoomedPlot(myLattice, spotDictionary, unassembledDataFile, processingFolder, runNumber, imageNumber, detectorDistance, pixelSize):
    goodPredictions_x = []
    goodPredictions_y = []
    goodPredictions_labels = []
    badPredictions_x = []
    badPredictions_y = []
    badPredictions_labels = []
    zoomedSpots_x = []
    zoomedSpots_y = []
    zoomedSpots_labels = []
                                               
    for spot in myLattice.refinedPredictedPattern:
        h = spot[0]
        k = spot[1]
        label = r'%d, %d $\rightarrow$ %d, %d'%(h, k, k, h)
        xPredicted = 874 + myLattice.refinedCenterXs[-1] + spot[10] * numpy.cos(spot[8]) 
        yPredicted = 874 + myLattice.refinedCenterYs[-1] + spot[10] * numpy.sin(spot[8]) 
        integrationFlag = 0
        for mySpotKey, mySpot in spotDictionary.items():  
            if mySpot.h == h and mySpot.k == k:
                try:
                    I = mySpot.integratedIntensity
                    if not numpy.isnan(I):
                        integrationFlag = 1
                except:
                    integrationFlag = 0
        if integrationFlag == 0:
            badPredictions_x.append(xPredicted)
            badPredictions_y.append(yPredicted)
            badPredictions_labels.append(label)
        else:
            goodPredictions_x.append(xPredicted)
            goodPredictions_y.append(yPredicted)
            goodPredictions_labels.append(label)
        if h == -13 and k == 11 or h == -11 and k == 13 or h == -13 and k == 2 or h == -11 and k == -2:
            if h == -11 and k == -2:                 # -11, -2 -> p -> -2, -11
                label = r'{(-2, -11)}'
            if h == -13 and k == 2:                  # 2, 11   -> p -> 11, 2
                label = r'{(11, 2)}'
            if h == -11 and k == 13:                 # -2, -11 -> p -> -11, -2
                label = r'{(-11, -2)}'
            if h == -13 and k == 11:                 # 11, 2   -> p -> 2, 11
                label = r'{(2, 11)}'
            zoomedSpots_x.append(xPredicted)
            zoomedSpots_y.append(yPredicted)
            zoomedSpots_labels.append(label)
            

    ### EXTRACT ASSEMBLED DATA FOR FINAL PLOTTING ###
    assembledData   = unassembledDataFile['/data/assembleddata0'] 
    assembledData   = numpy.asarray(assembledData, dtype=numpy.float32)         #### !!!!! ####
    
    ### INDEXED PREDICTED PATTERN FIGURE ###
    matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
    matplotlib.pyplot.title('%s'%myLattice.fileName, y=1.05)
    matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    
    matplotlib.pyplot.scatter(badPredictions_x, badPredictions_y, color='b', marker="o", linewidth='2', facecolors='none', s=600)
    for label, x, y in zip(badPredictions_labels, badPredictions_x, badPredictions_y):
        matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                   textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                   bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    matplotlib.pyplot.scatter(goodPredictions_x, goodPredictions_y,   color='r', marker="o", linewidth='2', facecolors='none', s=600)
    for label, x, y in zip(goodPredictions_labels, goodPredictions_x, goodPredictions_y):
        matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                   textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                   bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    
    matplotlib.pyplot.savefig("%s/r%s_Image_%s_Lattice_%s.png"%(processingFolder, runNumber, imageNumber, myLattice.latticeNumberInImage), fontsize = 20, dpi=96*4)
    matplotlib.pyplot.close()
    
    ### ZOOMED-IN FIGURE ###                
    matplotlib.pyplot.figure(figsize=(40,50), dpi=4*96, facecolor='w',frameon=True)  #width, height
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    
    ax1 = matplotlib.pyplot.subplot2grid((5,4), (0,0), rowspan=4, colspan=4)
    ax1.imshow(assembledData[250:1300, 700:1750], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
                 
    badPredictions_x = numpy.asarray(badPredictions_x)
    badPredictions_x = badPredictions_x - 700
    goodPredictions_x = numpy.asarray(goodPredictions_x)
    goodPredictions_x = goodPredictions_x - 700     
    
    badPredictions_y = numpy.asarray(badPredictions_y)
    badPredictions_y = badPredictions_y - 250
    goodPredictions_y = numpy.asarray(goodPredictions_y)
    goodPredictions_y = goodPredictions_y - 250
    
    ax1.scatter(badPredictions_x,  badPredictions_y,  color='b', marker="o", linewidth='2', facecolors='none', s=600)
    ax1.scatter(goodPredictions_x, goodPredictions_y, color='r', marker="o", linewidth='2', facecolors='none', s=600)
    ax1.scatter(874-700, 874-250, s=3600, marker='+', color='m')
    resolutionCircle = 7.0 #A
    myRadius = detectorDistance/pixelSize*math.tan(2*math.asin(myLattice.wavelength/(2*resolutionCircle)))
    circle = matplotlib.pyplot.Circle((874-700, 874-250), myRadius, linestyle='dashed', linewidth=2.0, color='b',fill=False)
    ax1.add_artist(circle)
    
    fs = 57
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.set_xlim([0, 1050])
    ax1.set_ylim([0, 1050])
    t1 = ax1.set_title('(a)', fontsize=fs, loc='left')                
    t1.set_position([0, 1.01])
    
    ax2 = matplotlib.pyplot.subplot2grid((5,4), (4,0), colspan=1)
    ax2.imshow(assembledData[zoomedSpots_y[0]-21:zoomedSpots_y[0]+22, zoomedSpots_x[0]-21:zoomedSpots_x[0]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, linestyle='dashed', color='b',fill=False)
    ax2.add_artist(circle)
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])
    t2_l = ax2.set_title('(b)', fontsize=fs, loc='left')                
    t2_r = ax2.set_title("%s"%(zoomedSpots_labels[0]), fontsize=fs)
    t2_l.set_position([0, 1.03])
    t2_r.set_position([0.5, 1.03])
    
    ax3 = matplotlib.pyplot.subplot2grid((5,4), (4,1), colspan=1)
    ax3.imshow(assembledData[zoomedSpots_y[1]-21:zoomedSpots_y[1]+22, zoomedSpots_x[1]-21:zoomedSpots_x[1]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, linestyle='dashed', color='b',fill=False)
    ax3.add_artist(circle)
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    t3_l = ax3.set_title('(c)', fontsize=fs, loc='left')
    t3_r = ax3.set_title("%s"%(zoomedSpots_labels[1]), fontsize=fs)
    t3_l.set_position([0, 1.03])
    t3_r.set_position([0.5, 1.03])
    
    ax4 = matplotlib.pyplot.subplot2grid((5,4), (4,2), colspan=1)
    ax4.imshow(assembledData[zoomedSpots_y[2]-21:zoomedSpots_y[2]+22, zoomedSpots_x[2]-21:zoomedSpots_x[2]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, linestyle='dashed', color='b',fill=False)
    ax4.add_artist(circle)
    ax4.set_xticklabels([])
    ax4.set_yticklabels([])
    t4_l = ax4.set_title('(d)', fontsize=fs, loc='left')
    t4_r = ax4.set_title("%s"%(zoomedSpots_labels[2]), fontsize=fs)
    t4_l.set_position([0, 1.03])
    t4_r.set_position([0.5, 1.03])
    
    ax5 = matplotlib.pyplot.subplot2grid((5,4), (4,3), colspan=1)
    ax5.imshow(assembledData[zoomedSpots_y[3]-21:zoomedSpots_y[3]+22, zoomedSpots_x[3]-21:zoomedSpots_x[3]+22], origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    circle = matplotlib.pyplot.Circle((21,21), 10, linewidth=1.0, linestyle='dashed', color='b',fill=False)
    ax5.add_artist(circle)
    ax5.set_xticklabels([])
    ax5.set_yticklabels([])
    t5_l = ax5.set_title('(e)', fontsize=fs, loc='left')
    t5_r = ax5.set_title("%s"%(zoomedSpots_labels[3]), fontsize=fs)
    t5_l.set_position([0, 1.03])
    t5_r.set_position([0.5, 1.03])
               
    matplotlib.pyplot.savefig("%s/r%s_Image_%s_Lattice_%s_zoomed.png"%(processingFolder, runNumber, imageNumber, myLattice.latticeNumberInImage), fontsize = 20, dpi=96*4)
    matplotlib.pyplot.savefig("%s/r%s_Image_%s_Lattice_%s_zoomed.pdf"%(processingFolder, runNumber, imageNumber, myLattice.latticeNumberInImage), fontsize = 20, dpi=96*4)
    matplotlib.pyplot.close()
