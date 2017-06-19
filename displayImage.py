# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 15:49:56 2016
@author: casadei_c
MEMBER OF diffractionImage CLASS
DISPLAY ***ASSEMBLED*** IMAGE ASSOCIATED TO EACH diffractionImage OBJECT.
"""
import matplotlib.pyplot
import h5py
import warnings

myDPI = 96
def displayImageFunction(self, folderName):
    warnings.filterwarnings("ignore")
    print 'Displaying image n %s: %s'%(self.imageNumber, self.fileName)
    matplotlib.pyplot.close()
    imgFile = h5py.File('%s/%s' %(folderName, self.fileName), 'r')
    assembledData = imgFile['data/assembleddata0']
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10, 10), dpi=myDPI, facecolor='w',frameon=True)
    matplotlib.pyplot.title('%s\nRun number: %s - Image number %s'%(self.fileName, self.runNumber, self.imageNumber), y=1.05)
    myAxesImageObject = matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 40, vmax = 100)
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.show()