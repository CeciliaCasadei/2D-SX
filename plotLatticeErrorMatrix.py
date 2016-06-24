# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 19:19:34 2016
@author: casadei_c
- plotLatticeErrorMatrixFunction:
  MEMBER OF Lattice CLASS
  CALLED BY orientationAndCellRefinement.py
  PLOT ATTRIBUTE latticeErrorMatrix
"""
import warnings
import matplotlib.pyplot

def plotLatticeErrorMatrixFunction(self, filename):    
    matplotlib.pyplot.close()        
    warnings.filterwarnings("ignore")
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10,10), dpi=80, facecolor='w',frameon=True)
    matplotlib.pyplot.title('%s\nImage number %s \nLattice: %d'%(self.fileName, self.imageNumber, self.latticeNumberInImage), fontsize=16, y=1.02)
    myAxesImageObject = matplotlib.pyplot.imshow(self.latticeErrorMatrix, origin='lower', interpolation='nearest')   
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    matplotlib.pyplot.gca().set_xlabel('Cell size', fontsize = 18, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel('In-plane rotation', fontsize = 18, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s'%filename)
    matplotlib.pyplot.close()