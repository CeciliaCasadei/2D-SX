# -*- coding: utf-8 -*-
import matplotlib.pyplot
import warnings



def plotGeometry(geoX, geoY):
    warnings.filterwarnings("ignore")
    matplotlib.pyplot.close()
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10, 10), dpi=96, facecolor='w',frameon=True)
    matplotlib.pyplot.title('X geometry', y=1.05)
    myAxesImageObject = matplotlib.pyplot.imshow(geoX, origin='lower', interpolation='nearest')
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.savefig('./geometryX.png')
    matplotlib.pyplot.close()
    myFigureObject  = matplotlib.pyplot.figure(figsize=(10, 10), dpi=96, facecolor='w',frameon=True)
    matplotlib.pyplot.title('Y geometry', y=1.05)
    myAxesImageObject = matplotlib.pyplot.imshow(geoY, origin='lower', interpolation='nearest')
    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.savefig('./geometryY.png')
    matplotlib.pyplot.close()