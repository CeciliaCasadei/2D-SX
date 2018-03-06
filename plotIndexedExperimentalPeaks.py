# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 15:04:22 2016
@author: casadei_c
MEMBER OF diffractionImage CLASS.
PLOT INDEXING RESULTS.
PLOT EXPERIMENTAL PEAKS (FROM peaks.txt)
AND INDEXED, CALCULATED PATTERN(S).
"""
import warnings
import os
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import numpy
import math
import pickle

def plotIndexedExperimentalPeaksFunction(self, 
                                         detectorDistance, 
                                         pixelSize, 
                                         resolutionRadii):
    outFolder = './Output_r%s/LatticeIndexing'%self.runNumber
    latticesFile = '%s/latticeDictionary_r%s_image_%s.pkl'%(outFolder, 
                                                            self.runNumber, 
                                                            self.imageNumber) 
    if os.path.exists(latticesFile):
        openPkl = open(latticesFile, 'rb')
        myLatticesDictionary = pickle.load(openPkl)
        openPkl.close()        
    else:
        return
        
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.close()
    myDPI = 96 
        
    x_data = self.orderedPeaksMatrix[:, 0]
    y_data = self.orderedPeaksMatrix[:, 1]
    
    myFigure = matplotlib.pyplot.figure(figsize=(25, 25), 
                                        dpi=myDPI, 
                                        facecolor = 'w')   
    myAxes = myFigure.add_subplot(1,1,1) 
    
    myAxes.scatter(x_data, y_data, color="red", marker="+")  
           
    totalNLattices = 0
    for i,j in myLatticesDictionary.items(): # Loop on lattices
        totalNLattices = totalNLattices + 1
    if totalNLattices == 1:
        myTextSize = 20
    elif totalNLattices == 2:
        myTextSize = 11
    elif totalNLattices == 3:
        myTextSize = 8
    else:
        myTextSize = 6        
        
    colorsList = ['k', 'c', 'b', 'g', 'r', 'm']
    nScatter = 0
    for i,j in myLatticesDictionary.items():
        myColor = colorsList[nScatter]
        myX = []
        myY = []
        myLabels = []
        
        for myRow in myLatticesDictionary['%s'%i].indexedPeaksTable:
            h_idx = myRow[0]
            h_idx = int(h_idx)
            k_idx = myRow[1]
            k_idx = int(k_idx)
            hk = '%d, %d'%(h_idx, k_idx)
            predictedRadius = myRow[9]
            predictedAzimuth = myRow[10]
            predicted_x = predictedRadius * math.cos(predictedAzimuth)
            predicted_y = predictedRadius * math.sin(predictedAzimuth)
            myX.append(predicted_x)
            myY.append(predicted_y)
            myLabels.append(hk)
            
        # REMOVE DUPLICATES
        myLabels_cleaned = []
        myX_cleaned = []
        myY_cleaned = []
        for i in range(0, len(myLabels)):
            label = myLabels[i]
            if label not in myLabels_cleaned:
                myLabels_cleaned.append(label)
                myX_cleaned.append(myX[i])
                myY_cleaned.append(myY[i])

        myAxes.scatter(myX_cleaned, 
                       myY_cleaned, 
                       color=myColor, 
                       marker="o", 
                       facecolors='none', 
                       s=70)
        for label, x, y in zip(myLabels_cleaned, myX_cleaned, myY_cleaned):
            matplotlib.pyplot.annotate(label, 
                                       xy = (x, y), 
                                       xytext = (-3, 3), 
                                       size = myTextSize, 
                                       textcoords = 'offset points', 
                                       ha = 'right', 
                                       va = 'bottom', 
                                       bbox = dict(boxstyle = 'round,pad=0.3', 
                                                   fc = 'yellow', 
                                                   alpha = 0.3, 
                                                   ec='none'))
        nScatter = nScatter + 1
        
    myRadii = []
    for i in resolutionRadii:
        i = float(i)
        myRadius = detectorDistance/pixelSize*numpy.tan(2*numpy.arcsin(self.wavelength/(2*i)))
        circle = matplotlib.pyplot.Circle((0,0),
                                           myRadius,
                                           linewidth=0.5, 
                                           color='b',
                                           fill=False)
        myAxes.add_artist(circle)
        myRadii.append(myRadius)
    
    matplotlib.pyplot.axhline(y=0, xmin=-870, xmax=870, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-870, ymax=870, linewidth=0.5, color = 'b')
    myAxes.set_title("Image %s: %s\nResolution circles: %.1f A, %.1f A, %.1f A\nNumber of lattices found: %d"%(self.imageNumber, 
                                                                                                               self.fileName, 
                                                                                                               resolutionRadii[0], 
                                                                                                               resolutionRadii[1], 
                                                                                                               resolutionRadii[2], 
                                                                                                               totalNLattices), 
                                                                                                               y=1.02, 
                                                                                                               fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    
    myAxes.set_xlim([-max(myRadii)-20,+max(myRadii)+20])
    myAxes.set_ylim([-max(myRadii)-20,+max(myRadii)+20])
    
    myAxes.set_xlabel("Detector pxls, x", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("Detector pxls, y", fontsize = 22, rotation = 'vertical')
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    
#    print 'Plotting figure n: %s (%s) \t %d lattice(s) found.'%(self.imageNumber, 
#                                                                self.fileName, 
#                                                                nScatter)
                                                                
    figuresPath = './Output_r%s/LatticeIndexing/Figures'%self.runNumber
    myFigure.savefig("%s/IndexedExperimentalPeaks_r%s_Image_%s.png"%(figuresPath, 
                                                                     self.runNumber, 
                                                                     self.imageNumber), 
                                                                     fontsize = 20, 
                                                                     dpi=myDPI*2)
    matplotlib.pyplot.close()