# -*- coding: utf-8 -*-
import matplotlib.pyplot
import warnings
import h5py
import numpy
import os


def plotGeometry(geoX, geoY):
    warnings.filterwarnings("ignore")
    outDirectory = './GeometryPlots'
    if not os.path.exists(outDirectory):
        os.mkdir(outDirectory)
    matplotlib.pyplot.close()
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10, 10), dpi=96, facecolor='w',frameon=True)
    matplotlib.pyplot.title('X geometry', y=1.05)
    myAxesImageObject = matplotlib.pyplot.imshow(geoX, origin='lower', interpolation='nearest')
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.savefig('%s/geometryX.png'%outDirectory)
    matplotlib.pyplot.close()
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10, 10), dpi=96, facecolor='w',frameon=True)
    matplotlib.pyplot.title('Y geometry', y=1.05)
    myAxesImageObject = matplotlib.pyplot.imshow(geoY, origin='lower', interpolation='nearest')
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.savefig('%s/geometryY.png'%outDirectory)
    matplotlib.pyplot.close()
    
    
def getGeometry():
    geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5' # same for all runs
    ### EXTRACT GEOMETRY ###
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    plotGeometry(xGeometry_np, yGeometry_np)
    
if __name__ == "__main__":
    print "\n**** CALLING plotGeometry ****"
    getGeometry()