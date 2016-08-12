# -*- coding: utf-8 -*
import numpy
import joblib
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

IsAtZero = joblib.load('./Output_runMerging/intensitiesAtZero/intensitiesAtZero.jbl')
print IsAtZero

a = 62.45
directCell = a * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2* numpy.pi * directCell.I     
reciprocalCellColumns = reciprocalCellRows.T

aStar_x = reciprocalCellColumns[0, 0]
aStar_y = reciprocalCellColumns[1, 0]
bStar_x = reciprocalCellColumns[0, 1]
bStar_y = reciprocalCellColumns[1, 1]

print aStar_x
print aStar_y
print bStar_x
print bStar_y



N = 17

n_s = range(0, 16+1)  # 17 values, 0 to 16
m_s = range(0, 16+1)  # 17 values, 0 to 16

h_s = range(-8, 8+1)  # 17 values, -8 to 8
k_s = range(-8, 8+1)  # 17 values, -8 to 8

xCart_s = []
yCart_s = []
patterson = []

hs_plot = []
ks_plot = []

q_xs = []
q_ys = []

for n in n_s:
    print n
    u = float(n)/N
    for m in m_s:
        v = float(m)/N
        
        pattersonSum = 0        
        for h in h_s:
            for k in k_s:
                nMatches = 0
                for braggRod in IsAtZero:
                    hRod = braggRod[0]
                    kRod = braggRod[1]
                    if (h == hRod and k == kRod) or (h == -hRod-kRod and k == hRod) or (h == kRod and k == -hRod-kRod) or (h == -hRod and k == -kRod) or (h == hRod+kRod and k == -hRod) or (h == -kRod and k == hRod+kRod):
                        pattersonSum = pattersonSum + braggRod[2] * numpy.cos(2*numpy.pi*(h*u + k*v))
                        nMatches = nMatches + 1
                        
                        ### PLOT ###
                        hs_plot.append(h)
                        ks_plot.append(k)
                        reciprocalVector = [h, k]*reciprocalCellRows
                        q_x = reciprocalVector[0,0]         # A^(-1)
                        q_y = reciprocalVector[0,1]         # A^(-1)    
                        
                        q_xs.append(q_x)
                        q_ys.append(q_y)
                        ###
                        if nMatches > 1:
                            print 'nMatches > 1'
                            
        xCart = u * a + v * a * numpy.cos(2*numpy.pi/3)           
        yCart = v * a * numpy.sin(2*numpy.pi/3)
        xCart_s.append(xCart)
        yCart_s.append(yCart)
        patterson.append(pattersonSum)
                  
        
fig = matplotlib.pyplot.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(xCart_s, yCart_s, patterson, cmap=cm.seismic, linewidth=0.2)
ax.set_zlim(-4500, 4500)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
matplotlib.pyplot.axis('off')
matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
matplotlib.pyplot.show()

myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
myAxes = myFigure.add_subplot(3,1,3)    # Axes object  -  one row, one column, first plot    
myAxes.scatter(xCart_s, yCart_s, color="red", marker=".")
myAxes.tick_params(axis='x', labelsize=14)
myAxes.tick_params(axis='y', labelsize=14)
myAxes.set_xlabel("x", fontsize = 22, rotation = 'horizontal')
myAxes.set_ylabel("y", fontsize = 22, rotation = 'vertical')
myAxes.axis('equal')

myAxes_2 = myFigure.add_subplot(3,1,1)    # Axes object  -  one row, one column, first plot    
myAxes_2.scatter(hs_plot, ks_plot, color="red", marker=".")
myAxes_2.tick_params(axis='x', labelsize=14)
myAxes_2.tick_params(axis='y', labelsize=14)
myAxes_2.set_xlabel("h", fontsize = 22, rotation = 'horizontal')
myAxes_2.set_ylabel("k", fontsize = 22, rotation = 'vertical')
myAxes_2.axis('equal')

myAxes_3 = myFigure.add_subplot(3,1,2)    # Axes object  -  one row, one column, first plot    
myAxes_3.scatter(q_xs, q_ys, color="red", marker=".")
myAxes_3.tick_params(axis='x', labelsize=14)
myAxes_3.tick_params(axis='y', labelsize=14)
myAxes_3.set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 22, rotation = 'horizontal')
myAxes_3.set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 22, rotation = 'vertical')
myAxes_3.axis('equal')

matplotlib.pyplot.savefig('./Output_runMerging/intensitiesAtZero/DFT_realReciprocalSpace.png')

xCart_s_extended = []
yCart_s_extended = []
patterson_extended = []

for i in range(0, len(xCart_s)):
    xCart_s_extended.append(xCart_s[i])
    yCart_s_extended.append(yCart_s[i])
    patterson_extended.append(patterson[i])
    
for i in range(0, len(xCart_s)):                # translate along real space vector a
    xCart_s_extended.append(xCart_s[i]+a)
    yCart_s_extended.append(yCart_s[i])
    patterson_extended.append(patterson[i])
    
for i in range(0, len(xCart_s)):                   # translate along real space vector b
    xCart_s_extended.append(xCart_s[i]+a*numpy.cos(2*numpy.pi/3))
    yCart_s_extended.append(yCart_s[i]+a*numpy.sin(2*numpy.pi/3))
    patterson_extended.append(patterson[i])
    
for i in range(0, len(xCart_s)):                   # translate along real space vector a and b
    xCart_s_extended.append(xCart_s[i]+a*numpy.cos(2*numpy.pi/3) + a)
    yCart_s_extended.append(yCart_s[i]+a*numpy.sin(2*numpy.pi/3))
    patterson_extended.append(patterson[i])
        
fig = matplotlib.pyplot.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(xCart_s_extended, yCart_s_extended, patterson_extended, cmap=cm.seismic, linewidth=0.2)
ax.set_zlim(-4500, 4500)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])
matplotlib.pyplot.axis('off')
matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
matplotlib.pyplot.show()     