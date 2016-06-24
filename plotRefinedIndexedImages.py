# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:03:29 2016
@author: casadei_c
MEMBER OF dffractionImage CLASS
FOR EACH DIFFRACTION IMAGE
PLOT EXPERIMENTAL PEAKS (FROM peaks.txt) AND REFINED, INDEXED PREDICTED LATTICE(S)
"""
import warnings
import matplotlib.pyplot
import pickle
import math
import os

def plotRefinedIndexedImagesFunction(self, resolutionRadii):
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.close()
    myDPI = 96 
        
    x_data = self.orderedPeaksMatrix[:, 0]
    y_data = self.orderedPeaksMatrix[:, 1]
    
    myFigure = matplotlib.pyplot.figure(figsize=(25, 25), dpi=myDPI, facecolor = 'w')   # Figure object
    myAxes = myFigure.add_subplot(1,1,1) # Axes object
    
    myAxes.scatter(x_data, y_data, color="red", marker="+")  
    
    colorsList = ['k', 'c', 'b', 'g', 'r', 'm']
    fLattices = open('./Output_r%s/OrientationAndCellRefinement/r%s_refinedLatticesDictionary.pkl'%(self.runNumber, self.runNumber), 'rb')
    myDataLattices = pickle.load(fLattices)
    fLattices.close()
    totalNLattices = 0
    for i,j in myDataLattices.items(): # Loop on lattices
        if myDataLattices['%s'%i].imageNumber == self.imageNumber: # Lattice belongs to this image
            totalNLattices = totalNLattices + 1
    if totalNLattices == 1:
        myTextSize = 20
    elif totalNLattices == 2:
        myTextSize = 11
    elif totalNLattices == 3:
        myTextSize = 8
    else:
        myTextSize = 6      
            
    nLattices = 0
    for i,j in myDataLattices.items(): # Loop on lattices
        if myDataLattices['%s'%i].imageNumber == self.imageNumber: # Lattice belongs to this image
            myColor = colorsList[nLattices]
            myX = []
            myY = []
            myLabels = []
            nLattices = nLattices + 1
            detectorDistance = myDataLattices['%s'%i].detectorDistance
            pixelSize = myDataLattices['%s'%i].pixelSize
            refinedPredictedPattern = myDataLattices['%s'%i].refinedPredictedPattern
            for myRow in refinedPredictedPattern:
                h_idx = myRow[0]
                h_idx = int(h_idx)
                k_idx = myRow[1]
                k_idx = int(k_idx)
                hk = '%d, %d'%(h_idx, k_idx)
                predictedRadius  = myRow[10]
                predictedAzimuth = myRow[8]
                predicted_x = predictedRadius * math.cos(predictedAzimuth)
                predicted_y = predictedRadius * math.sin(predictedAzimuth)
                for indexedPeaks in myDataLattices['%s'%i].indexedPeaksTable:
                    if indexedPeaks[0] == h_idx and indexedPeaks[1] == k_idx:
                        myX.append(predicted_x)
                        myY.append(predicted_y)
                        myLabels.append(hk)
            myAxes.scatter(myX, myY, color=myColor, marker="o", facecolors='none', s=70)
            for label, x, y in zip(myLabels, myX, myY):
                matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = myTextSize, 
                                           textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                           bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    
    if nLattices == 0:
        return
                                                                                 
    myRadii = []
    for i in resolutionRadii:
        i = float(i)
        myRadius = detectorDistance/pixelSize*math.tan(2*math.asin(self.wavelength/(2*i)))
        circle = matplotlib.pyplot.Circle((0,0),myRadius,linewidth=0.5, color='b',fill=False)
        myAxes.add_artist(circle)
        myRadii.append(myRadius)
    
    matplotlib.pyplot.axhline(y=0, xmin=-870, xmax=870, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-870, ymax=870, linewidth=0.5, color = 'b')
    myAxes.set_title("Image %s: %s\nResolution circles: %.1f A, %.1f A, %.1f A\nNumber of lattices found: %d"%(self.imageNumber, self.fileName, resolutionRadii[0], resolutionRadii[1], resolutionRadii[2], totalNLattices), y=1.02, fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    
    myAxes.set_xlim([-max(myRadii)-20,+max(myRadii)+20])
    myAxes.set_ylim([-max(myRadii)-20,+max(myRadii)+20])
    
    myAxes.set_xlabel("Detector pxls, x", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("Detector pxls, y", fontsize = 22, rotation = 'vertical')
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    
    if not os.path.exists('./Output_r%s/OrientationAndCellRefinement/Figures'%self.runNumber):
        os.mkdir('./Output_r%s/OrientationAndCellRefinement/Figures'%self.runNumber)   
    myFigure.savefig("./Output_r%s/OrientationAndCellRefinement/Figures/IndexedExperimentalPeaks_Refined_r%s_Image_%s.png"%(self.runNumber, self.runNumber, self.imageNumber), fontsize = 20, dpi=myDPI*2)
    matplotlib.pyplot.close()
    print 'Plotting image %s'%self.fileName