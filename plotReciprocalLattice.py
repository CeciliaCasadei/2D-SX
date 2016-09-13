# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:43:29 2016
@author: casadei_c
PRODUCE RECIPROCAL LATTICE FIGURE.
"""
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import pickle
import warnings 

def plotReciprocalLatticeFunction(cellSize, outFolder):
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.ion()
    fRead = open('%s/ReferenceReciprocalLattice/reciprocalLattice_cellSize_%.3f.pkl'%(outFolder, cellSize), 'rb')
    myData = pickle.load(fRead)
    fRead.close()
    
    x_data = []
    y_data = []
    myLabels = []
    for i in myData:
        x_data.append(i[2])
        y_data.append(i[3])        
        h_idx = i[0]
        h_idx = int(h_idx)
        k_idx = i[1]
        k_idx = int(k_idx)
        hk = '%d, %d'%(h_idx, k_idx)
        myLabels.append(hk)
    
    myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
    myAxes = myFigure.add_subplot(1,1,1)    # Axes object  -  one row, one column, first plot    
    myAxes.scatter(x_data, y_data, color="red", marker=".")
    for label, x, y in zip(myLabels, x_data, y_data):
        matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 9, 
                                   textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                   bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    
    myAxes.set_title("Reciprocal lattice, reference orientation and cell size", y=1.05, fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    
    matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-1, ymax=1, linewidth=0.5, color = 'b')
    
    myAxes.set_xlim([min(x_data)-0.1,max(x_data)+0.1])
    myAxes.set_ylim([min(y_data)-0.1,max(y_data)+0.1])
    
    myAxes.set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 22, rotation = 'vertical')
    
    matplotlib.pyplot.show()
    myFigure.savefig("%s/ReferenceReciprocalLattice/reciprocalLattice_cellSize_%.3f.png"%(outFolder, cellSize), fontsize = 20)
    print 'Reciprocal lattice plot saved in %s/ReferenceReciprocalLattice/reciprocalLattice_cellSize_%.3f.png'%(outFolder, cellSize)