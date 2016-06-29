# -*- coding: utf-8 -*-
class diffractionSpot:
    
    def __init__(self, n, h, k, 
                 boxMatrix, maskBoxMatrix, bgSubtractedBoxMatrix, peakDetectionMask, xLeft, yDown):
        self.n = n
        self.h = h
        self.k = k
        self.boxMatrix = boxMatrix
        self.maskBoxMatrix = maskBoxMatrix
        self.bgSubtractedBoxMatrix = bgSubtractedBoxMatrix
        self.peakDetectionMask = peakDetectionMask
        self.xLeft = xLeft
        self.yDown = yDown        
        # detectedSpot objects objects are generated in unassembledImageProcessing


        
    def setDistanceFromPrediction(self, distanceFromPrediction):
        self.distanceFromPrediction = distanceFromPrediction
        
        
        
    def setBoxIndices(self, i, j):
        self.i = i
        self.j = j
        
        
        
    def setFinalBoxIndices(self, i, j):
        self.iFinal = i
        self.jFinal = j
        
        
        
    def integrateSpot(self):
        import integrateSpot
        integratedIntensity = integrateSpot.integrate(self)
        self.integratedIntensity = integratedIntensity
        
        
        
    def spotProcessingPlot(self, runNumber, imageNumber, latticeNumberInImage, boxWidth):
        import warnings
        import matplotlib.pyplot
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        
        warnings.filterwarnings("ignore")  
        matplotlib.pyplot.ioff                                                                 
        myFigure = matplotlib.pyplot.figure()
        myFigure.suptitle('Peak %d \nh = %d    k = %d'%(self.n, self.h, self.k), fontsize=12)
        
        ax1 = myFigure.add_subplot(321)
        ax1.tick_params(axis='both', which='major', labelsize=6)
        ax1.set_title('UNASSEMBLED RAW DATA SECTOR', fontsize=5)
        im1 = ax1.imshow(self.boxMatrix, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)     
        divider1 = make_axes_locatable(ax1)    
        cax1 = divider1.append_axes("right", size="10%", pad=0.05)  
        cb1 = matplotlib.pyplot.colorbar(im1, cax=cax1)
        for t in cb1.ax.get_yticklabels():
            t.set_fontsize(6)    
    
        ax2 = myFigure.add_subplot(322)
        ax2.tick_params(axis='both', which='major', labelsize=6)
        ax2.set_title('MODULE MASK', fontsize=5)
        im2 = ax2.imshow(self.maskBoxMatrix, origin='lower', interpolation='nearest', vmin = 0, vmax = 1) 
        im2.set_cmap('Blues')   
        
        ax3 = myFigure.add_subplot(323)   
        ax3.tick_params(axis='both', which='major', labelsize=6)             
        ax3.set_title('BACKGROUND SUBTRACTED', fontsize=5)
        im3 = ax3.imshow(self.bgSubtractedBoxMatrix, origin='lower', interpolation='nearest', vmin = 0, vmax = 100) 
        divider3 = make_axes_locatable(ax3)    
        cax3 = divider3.append_axes("right", size="10%", pad=0.05)  
        cb3 = matplotlib.pyplot.colorbar(im3, cax=cax3) 
        for t in cb3.ax.get_yticklabels():
            t.set_fontsize(6)    
        
        ax4 = myFigure.add_subplot(324)
        ax4.tick_params(axis='both', which='major', labelsize=6)
        ax4.set_title('CONNECTED PEAKS', fontsize=5)
        im4 = ax4.imshow(self.peakDetectionMask, origin='lower', interpolation='nearest', vmin = 0, vmax = 1) 
        matplotlib.pyplot.axhline(y=self.i, xmin=0, xmax=boxWidth, linewidth=0.1, color = 'r')
        matplotlib.pyplot.axvline(x=self.j, ymin=0, ymax=boxWidth, linewidth=0.1, color = 'r')
        if self.iFinal >= 0 and self.jFinal >= 0:
            matplotlib.pyplot.axhline(y=self.iFinal, xmin=0, xmax=boxWidth, linewidth=0.1, color = 'b')
            matplotlib.pyplot.axvline(x=self.jFinal, ymin=0, ymax=boxWidth, linewidth=0.1, color = 'b')
        im4.set_cmap('Blues')

        ax5 = myFigure.add_subplot(325)   
        ax5.tick_params(axis='both', which='major', labelsize=6)             
        ax5.set_title('INTEGRATION CIRCLE', fontsize=5)
        im5 = ax5.imshow(self.bgSubtractedBoxMatrix, origin='lower', interpolation='nearest', vmin = 0, vmax = 100) 
        if self.iFinal >= 0 and self.jFinal >= 0:
            matplotlib.pyplot.axhline(y=self.iFinal, xmin=0, xmax=boxWidth, linewidth=0.5, color = 'w')
            matplotlib.pyplot.axvline(x=self.jFinal, ymin=0, ymax=boxWidth, linewidth=0.5, color = 'w')
            circle = matplotlib.pyplot.Circle((self.jFinal, self.iFinal), self.integrationRadius, linewidth=0.5, color='w', fill = False)
            ax5.add_artist(circle)
        divider5 = make_axes_locatable(ax5)    
        cax5 = divider5.append_axes("right", size="10%", pad=0.05)  
        cb5 = matplotlib.pyplot.colorbar(im5, cax=cax5)         
        for t in cb5.ax.get_yticklabels():
            t.set_fontsize(6)  

        if hasattr(self, 'integrationMask'):
            ax6 = myFigure.add_subplot(326)
            ax6.tick_params(axis='both', which='major', labelsize=6)
            ax6.set_title('INTEGRATION MASK (FLAG %d)'%self.integrationFlag, fontsize=5)
            im6 = ax6.imshow(self.integrationMask, origin='lower', interpolation='nearest', vmin = 0, vmax = 1) 
            matplotlib.pyplot.axhline(y=self.i, xmin=0, xmax=boxWidth, linewidth=0.1, color = 'r')
            matplotlib.pyplot.axvline(x=self.j, ymin=0, ymax=boxWidth, linewidth=0.1, color = 'r')
            if self.iFinal >= 0 and self.jFinal >= 0:
                matplotlib.pyplot.axhline(y=self.iFinal, xmin=0, xmax=boxWidth, linewidth=0.1, color = 'b')
                matplotlib.pyplot.axvline(x=self.jFinal, ymin=0, ymax=boxWidth, linewidth=0.1, color = 'b')                
            im6.set_cmap('BuPu')
     
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.subplots_adjust(top=0.9)      
        matplotlib.pyplot.savefig('./Output_r%s/UnassembledImageProcessing/BgSubtractionAndPeakSearchPlots/r%s_img%s_lattice%s_spot%d.png'%(runNumber, runNumber, imageNumber, latticeNumberInImage, self.n), dpi = 4*96)
        matplotlib.pyplot.close()        