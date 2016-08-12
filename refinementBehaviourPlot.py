# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 16:55:35 2016
@author: casadei_c
MEMBER OF Lattice CLASS
PLOT BEHAVIOUR OF IN-PLANE ORIENTATION, CELL SIZE, CENTER COORDINATES, 
NUMBER OF DETECTED SPOTS AND LATTICE ERROR 
(OBSERVED TO PREDICTED SPOTS DISCREPANCY)
DURING REFINEMENT ITERATIONS IN IMAGE PROCESSING.
"""
import warnings
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import numpy

def refinementBehaviourPlotFunction(self, folderName, myString):
    try:
        nIterations = len(self.refinedLatticeOrientations) - 1
        xRange = range(0, nIterations+1)
        
        warnings.filterwarnings("ignore")                                                                   
        myFigure, (ax1, ax2, ax3, ax4, ax5, ax6) = matplotlib.pyplot.subplots(6, 1)    
        myFigure.suptitle('Run: %s - Image: %s - Lattice n: %s'%(self.runNumber, self.imageNumber, self.latticeNumberInImage), fontsize=15)
        
        ax1.set_title('ORIENTATION REFINEMENT', fontsize=13)
        ax1.plot(xRange, self.refinedLatticeOrientations, marker='o', linestyle='-', color='r')
        ax1.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.refinedLatticeOrientations)-min(self.refinedLatticeOrientations))/2
        if yStep < 0.000001:
            ax1.yaxis.set_ticks(numpy.arange(min(self.refinedLatticeOrientations)-0.0001,max(self.refinedLatticeOrientations)+0.0002, 0.0001))
            ax1.set_ylim([min(self.refinedLatticeOrientations) - 0.0001, max(self.refinedLatticeOrientations) + 0.0001])
        else:
            ax1.yaxis.set_ticks(numpy.arange(min(self.refinedLatticeOrientations),max(self.refinedLatticeOrientations)+yStep,yStep))
            ax1.set_ylim([min(self.refinedLatticeOrientations) - yStep/2, max(self.refinedLatticeOrientations) + yStep/2])
    
        ax2.set_title('CELL REFINEMENT', fontsize=13)
        ax2.plot(xRange, self.refinedCellSizes, marker='o', linestyle='-', color='r')  
        ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.refinedCellSizes)-min(self.refinedCellSizes))/2
        if yStep < 0.000001:
            ax2.yaxis.set_ticks(numpy.arange(min(self.refinedCellSizes)-0.0001,max(self.refinedCellSizes)+0.0002, 0.0001))
            ax2.set_ylim([min(self.refinedCellSizes) - 0.0001, max(self.refinedCellSizes) + 0.0001])
        else:
            ax2.yaxis.set_ticks(numpy.arange(min(self.refinedCellSizes),max(self.refinedCellSizes)+yStep,yStep))
            ax2.set_ylim([min(self.refinedCellSizes) - yStep/2, max(self.refinedCellSizes) + yStep/2])
    
        ax3.set_title('X CENTER', fontsize=13)
        ax3.plot(xRange, self.refinedCenterXs, marker='o', linestyle='-', color='r')
        ax3.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.refinedCenterXs)-min(self.refinedCenterXs))/2
        if yStep != 0:
            ax3.yaxis.set_ticks(numpy.arange(min(self.refinedCenterXs),max(self.refinedCenterXs)+yStep,yStep))
            ax3.set_ylim([min(self.refinedCenterXs) - yStep/2, max(self.refinedCenterXs) + yStep/2])
        else:
            ax3.yaxis.set_ticks(numpy.arange(min(self.refinedCenterXs)-1,max(self.refinedCenterXs)+2, 1))
            ax3.set_ylim([min(self.refinedCenterXs) - 1, max(self.refinedCenterXs) + 1])
        
        ax4.set_title('Y CENTER', fontsize=13)
        ax4.plot(xRange, self.refinedCenterYs, marker='o', linestyle='-', color='r')
        ax4.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.refinedCenterYs)-min(self.refinedCenterYs))/2
        if yStep != 0:
            ax4.yaxis.set_ticks(numpy.arange(min(self.refinedCenterYs),max(self.refinedCenterYs)+yStep,yStep))
            ax4.set_ylim([min(self.refinedCenterYs) - yStep/2, max(self.refinedCenterYs) + yStep/2])
        else:
            ax4.yaxis.set_ticks(numpy.arange(min(self.refinedCenterYs)-1,max(self.refinedCenterYs)+2, 1))
            ax4.set_ylim([min(self.refinedCenterYs) - 1, max(self.refinedCenterYs) + 1])
        
        ax5.set_title('N DETECTED SPOTS', fontsize=13)
        ax5.plot(xRange, self.nDetectedAndMatchedPeaks, marker='o', linestyle='-', color='r')   
        ax5.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.nDetectedAndMatchedPeaks)-min(self.nDetectedAndMatchedPeaks))/2
        if yStep != 0:
            ax5.yaxis.set_ticks(numpy.arange(min(self.nDetectedAndMatchedPeaks),max(self.nDetectedAndMatchedPeaks)+yStep,yStep))
            ax5.set_ylim([min(self.nDetectedAndMatchedPeaks) - yStep/2, max(self.nDetectedAndMatchedPeaks) + yStep/2])
        
        ax6.set_title('AVERAGE (RELATIVE) ERROR PER PEAK', fontsize=13)
        ax6.plot(xRange, self.avgLatticeErrors, marker='o', linestyle='-', color='r')   
        ax6.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off') 
        yStep = (max(self.avgLatticeErrors)-min(self.avgLatticeErrors))/2
        ax6.yaxis.set_ticks(numpy.arange(min(self.avgLatticeErrors),max(self.avgLatticeErrors)+yStep,yStep))
        ax6.set_ylim([min(self.avgLatticeErrors) - yStep/2, max(self.avgLatticeErrors) + yStep/2])
         
        
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.subplots_adjust(top=0.9)
        
        matplotlib.pyplot.savefig('%s/r%s_img%s_lattice%s_refinementBehaviour%s.png'%(folderName, self.runNumber, self.imageNumber, self.latticeNumberInImage, myString), dpi = 96)
        matplotlib.pyplot.close()
        
    except:
        print 'EXCEPT'
        nIterations = len(self.refinedLatticeOrientations) - 1
        xRange = range(0, nIterations+1)
        
        warnings.filterwarnings("ignore")                                                                   
        myFigure, (ax1, ax2, ax3, ax4, ax5, ax6) = matplotlib.pyplot.subplots(6, 1)    
        myFigure.suptitle('Run: %s - Image: %s - Lattice n: %s'%(self.runNumber, self.imageNumber, self.latticeNumberInImage), fontsize=15)
        
        ax1.set_title('ORIENTATION REFINEMENT', fontsize=13)
        ax1.plot(xRange, self.refinedLatticeOrientations, marker='o', linestyle='-', color='r')

        ax2.set_title('CELL REFINEMENT', fontsize=13)
        ax2.plot(xRange, self.refinedCellSizes, marker='o', linestyle='-', color='r')  

        ax3.set_title('X CENTER', fontsize=13)
        ax3.plot(xRange, self.refinedCenterXs, marker='o', linestyle='-', color='r')

        ax4.set_title('Y CENTER', fontsize=13)
        ax4.plot(xRange, self.refinedCenterYs, marker='o', linestyle='-', color='r')

        ax5.set_title('N DETECTED SPOTS', fontsize=13)
        ax5.plot(xRange, self.nDetectedAndMatchedPeaks, marker='o', linestyle='-', color='r')   
   
        ax6.set_title('AVERAGE (RELATIVE) ERROR PER PEAK', fontsize=13)
        ax6.plot(xRange, self.avgLatticeErrors, marker='o', linestyle='-', color='r')   
     
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.subplots_adjust(top=0.9)
        
        matplotlib.pyplot.savefig('%s/EXCEPT-r%s_img%s_lattice%s_refinementBehaviour%s.png'%(folderName, self.runNumber, self.imageNumber, self.latticeNumberInImage, myString), dpi = 96)
        matplotlib.pyplot.close()