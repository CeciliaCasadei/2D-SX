# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 13:24:08 2016
@author: casadei_c
PLOT LOCAL BACKGROUND PLANE (ALSO INTERACTIVE MODE, ONLY STATIC IMAGE IS SAVED)
"""
import os
import matplotlib.pyplot
import warnings
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def plotBgPlane(myX, myY, myZ,
                myData,
                C,
                directoryName, runNumber, imageNumber, latticeNumberInImage, nPeak):
                    
    warnings.filterwarnings("ignore") 
    showFlag = 0
    if not os.path.exists('%s'%directoryName):
        os.mkdir('%s'%directoryName)
     
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    myAxes.scatter(myData[:,0], myData[:,1], myData[:,2], c='r', marker='o', s=6)
    mySurface = myAxes.plot_surface(myX, myY, myZ, rstride=1, cstride=1, cmap=cm.coolwarm,
                                    linewidth=0, antialiased=False, alpha = 0.3)
    myFigure.colorbar(mySurface, shrink=0.5, aspect=5)
  
    if C[1] >= 0:
        s1 = '+'
    else:
        s1 = '-'
    if C[2] >= 0:
        s2 = '+'
    else:
        s2 = '-'        
    matplotlib.pyplot.title('z = %.3fx %s %.3fy %s %.2f'%(C[0], s1, abs(C[1]), s2, abs(C[2])), fontsize = 10, y=1.08)
    
    myAxes.zaxis.set_major_locator(LinearLocator(10))
    myAxes.zaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    myAxes.set_zlim(0, 40)
    myAxes.set_xlabel('X')
    myAxes.set_ylabel('Y')
    myAxes.set_zlabel('Z')  
    myAxes.axis('equal')
    myAxes.axis('tight')
    myAxes.tick_params(labelsize=8)
   
    matplotlib.pyplot.savefig('%s/BgPlane_r%s_img%s_lattice%s_peak%s.png'
                               %(directoryName, runNumber, imageNumber, latticeNumberInImage, nPeak))
    if showFlag == 1:
        matplotlib.pyplot.show()
    matplotlib.pyplot.close()   