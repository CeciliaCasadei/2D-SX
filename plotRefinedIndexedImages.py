# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 10:03:29 2016
@author: casadei_c
MEMBER OF dffractionImage CLASS
FOR EACH DIFFRACTION IMAGE
PLOT EXPERIMENTAL PEAKS (FROM peaks.txt) AND REFINED, INDEXED PREDICTED LATTICE(S)
"""
import warnings
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import pickle
import math
import os
import numpy
import h5py

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
        
            myAxes.scatter(myX_cleaned, myY_cleaned, color=myColor, marker="o", facecolors='none', s=70)
            for label, x, y in zip(myLabels_cleaned, myX_cleaned, myY_cleaned):
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
    



def plotRefinedIndexedImagesFunction_imageOverlap(self, resolutionRadii):
    from matplotlib import rcParams
    rcParams['xtick.direction'] = 'in'
    rcParams['ytick.direction'] = 'in'
    if not os.path.exists('./Output_r%s/OrientationAndCellRefinement/Figures_imageOverlap'%self.runNumber):
        os.mkdir('./Output_r%s/OrientationAndCellRefinement/Figures_imageOverlap'%self.runNumber)   
                             
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.close()
    myDPI = 96 
    cellSize = 62.45
    resolutionLimit = 7.0
    center_x = 874
    center_y = 874
    
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I                                                   # A^(-1)
                                
    fLattices = open('./Output_r%s/OrientationAndCellRefinement/r%s_refinedLatticesDictionary.pkl'%(self.runNumber, self.runNumber), 'rb')
    myDataLattices = pickle.load(fLattices)
    fLattices.close()
    
    unassembledDataFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1/%s'%(self.runNumber, self.fileName), 'r')
    unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
    assembledData   = unassembledDataFile['/data/assembleddata0'] 
    
    unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)       #### !!!!! ####
    assembledData   = numpy.asarray(assembledData, dtype=numpy.float32)         #### !!!!! ####
                                                            
    # LOOP ON LATTICES
    for i,j in myDataLattices.items(): 
        if myDataLattices['%s'%i].imageNumber == self.imageNumber: # Lattice belongs to this image
            myLattice = myDataLattices['%s'%i]
            
            # PLOT ASSEMBLED IMAGE       
            matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
            matplotlib.pyplot.title('%s'%self.fileName, y=1.05)
            matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100, cmap='Greys')
    
            predicted_X = []
            predicted_Y = []
            matched_X = []
            matched_Y = []
                      
            detectorDistance = myLattice.detectorDistance
            pixelSize = myLattice.pixelSize
            refinedPredictedPattern = myLattice.refinedPredictedPattern
            latticeNumberInImage = myLattice.latticeNumberInImage
            
            # LOOP ON WHOLE PREDICTED PATTERN
            for myRow in refinedPredictedPattern: 
                
                h_idx = myRow[0]
                h_idx = int(h_idx)
                k_idx = myRow[1]
                k_idx = int(k_idx)
                
                reciprocalVector = [h_idx, k_idx]*reciprocalCellRows
                q_x = reciprocalVector[0,0]         # A^(-1)
                q_y = reciprocalVector[0,1]         # A^(-1)
                q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
                resolution = 2* numpy.pi / q        # A
                if resolution < resolutionLimit:
                    continue
                
                predictedRadius  = myRow[10]
                predictedAzimuth = myRow[8]
                predicted_x = center_x + predictedRadius * math.cos(predictedAzimuth)
                predicted_y = center_y + predictedRadius * math.sin(predictedAzimuth)
                
                # LOOP ON MATCHED PEAKS
                match_label = 0
                for indexedPeaks in myLattice.indexedPeaksTable:
                    if indexedPeaks[0] == h_idx and indexedPeaks[1] == k_idx:
                        match_label = 1
                        matched_X.append(predicted_x)
                        matched_Y.append(predicted_y)
                        break
                        
                if match_label == 0:
                    predicted_X.append(predicted_x)
                    predicted_Y.append(predicted_y)
            matplotlib.pyplot.scatter(predicted_X, predicted_Y, color='b', marker="o", linewidth='2', facecolors='none', s=600)
            matplotlib.pyplot.scatter(matched_X,   matched_Y,   color='r', marker="o", linewidth='2', facecolors='none', s=600)
                                                   
            myRadii = []
            for i in resolutionRadii:
                i = float(i)
                myRadius = detectorDistance/pixelSize*math.tan(2*math.asin(self.wavelength/(2*i)))
                circle = matplotlib.pyplot.Circle((center_x, center_y),myRadius,linewidth=1.0, color='b',fill=False)
                matplotlib.pyplot.gca().add_artist(circle)
                myRadii.append(myRadius)
    
            matplotlib.pyplot.axhline(y=center_y, xmin=-870, xmax=870, linewidth=1.0, color = 'b')
            matplotlib.pyplot.axvline(x=center_x, ymin=-870, ymax=870, linewidth=1.0, color = 'b')
            matplotlib.pyplot.gca().set_title("Image %s: %s\nResolution circles: %.1f A, %.1f A, %.1f A"%(self.imageNumber, self.fileName, resolutionRadii[0], resolutionRadii[1], resolutionRadii[2]), y=1.02, fontsize = 20)

            matplotlib.pyplot.gca().set_xticklabels([])
            matplotlib.pyplot.gca().set_yticklabels([])
            
            matplotlib.pyplot.gca().set_xlim([center_x-max(myRadii)-20, center_x+max(myRadii)+20])
            matplotlib.pyplot.gca().set_ylim([center_y-max(myRadii)-20, center_y+max(myRadii)+20])
            
            matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
                      
            matplotlib.pyplot.savefig("./Output_r%s/OrientationAndCellRefinement/Figures_imageOverlap/Overlap_IndexedExperimentalPeaks_Refined_r%s_Image_%s_Lattice_%d.png"%(self.runNumber, self.runNumber, self.imageNumber, latticeNumberInImage), fontsize = 20, dpi=myDPI*4)
            matplotlib.pyplot.close()
            print 'Image %s - N %s - Lattice %s - N matched %d'%(self.fileName, self.imageNumber, latticeNumberInImage, len(matched_X))