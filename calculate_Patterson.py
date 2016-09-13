# -*- coding: utf-8 -*
import numpy
import getopt
import sys
import joblib
import matplotlib
from matplotlib import cm
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def calculate_PattersonFunction(myArguments):
    inputFolder = './Output_runMergingVsModel'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder="])
    except getopt.GetoptError:
        print 'Usage: python merging.py --inputFolder <inputFolder>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python merging.py --inputFolder <inputFolder>'
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
            
            
    IsAtZero = joblib.load('%s/intensitiesAtZero/intensitiesAtZero.jbl'%inputFolder)
    print IsAtZero
    
    a = 62.45
    directCell = a * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I     
   
    N = 17
    
    n_s = range(0, 16+1)  # 17 values, 0 to 16
    m_s = range(0, 16+1)  # 17 values, 0 to 16
    
    h_s = range(-8, 8+1)  # 17 values, -8 to 8
    k_s = range(-8, 8+1)  # 17 values, -8 to 8
    
    xCart_s = []
    yCart_s = []
    patterson = []
        
    for n in n_s:
        print n
        u = float(n)/N
        for m in m_s:
            v = float(m)/N
            
            hs_plot = []
            ks_plot = []
    
            q_xs = []
            q_ys = []
            
            Is_plot= []
            
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
                            Is_plot.append(braggRod[2])
                            ###
                            
                            if nMatches > 1:
                                print 'nMatches > 1'
                                
            xCart = u * a + v * a * numpy.cos(2*numpy.pi/3)           
            yCart = v * a * numpy.sin(2*numpy.pi/3)
            xCart_s.append(xCart)
            yCart_s.append(yCart)
            patterson.append(pattersonSum)
            
            
    ############FIGURES###################
                
    # 3D plot of I(qx, qy) at qRod = 0        
    fig = matplotlib.pyplot.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(q_xs, q_ys, Is_plot, cmap=cm.seismic, linewidth=0.2)
    ax.set_zticklabels([])
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    matplotlib.pyplot.show()
                      
    # 3D plot of patterson function vs real space        
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
    
    # 2D plots of measured h, k; measured qx, qy and x, y point where patterson was calculated
    myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
    
    myAxes_1 = myFigure.add_subplot(3,1,1)    # Axes object  -  one row, one column, first plot    
    myAxes_1.scatter(hs_plot, ks_plot, color="red", marker=".")
    myAxes_1.tick_params(axis='x', labelsize=14)
    myAxes_1.tick_params(axis='y', labelsize=14)
    myAxes_1.set_xlabel("h", fontsize = 22, rotation = 'horizontal')
    myAxes_1.set_ylabel("k", fontsize = 22, rotation = 'vertical')
    myAxes_1.axis('equal')
    
    myAxes_2 = myFigure.add_subplot(3,1,2)    # Axes object  -  one row, one column, first plot    
    myAxes_2.scatter(q_xs, q_ys, color="red", marker=".")
    myAxes_2.tick_params(axis='x', labelsize=14)
    myAxes_2.tick_params(axis='y', labelsize=14)
    myAxes_2.set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 22, rotation = 'horizontal')
    myAxes_2.set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 22, rotation = 'vertical')
    myAxes_2.axis('equal')
    
    myAxes_3 = myFigure.add_subplot(3,1,3)    # Axes object  -  one row, one column, first plot    
    myAxes_3.scatter(xCart_s, yCart_s, color="red", marker=".")
    myAxes_3.tick_params(axis='x', labelsize=14)
    myAxes_3.tick_params(axis='y', labelsize=14)
    myAxes_3.set_xlabel("x", fontsize = 22, rotation = 'horizontal')
    myAxes_3.set_ylabel("y", fontsize = 22, rotation = 'vertical')
    myAxes_3.axis('equal')
        
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/real_reciprocalSpace.png'%inputFolder)
    matplotlib.pyplot.close()
    
    # 3D plot of translation-extended patterson.
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
    
    # Plot of measured qx, qy points, with labels.
    myLabels = []
    for i in range(0, len(hs_plot)):   
        h_idx = hs_plot[i]
        h_idx = int(h_idx)
        k_idx = ks_plot[i]
        k_idx = int(k_idx)
        hk = '%d, %d'%(h_idx, k_idx)
        myLabels.append(hk)
        
    myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
    myAxes = myFigure.add_subplot(1,1,1)    # Axes object  -  one row, one column, first plot    
    myAxes.scatter(q_xs, q_ys, color="red", marker=".")
    for label, x, y in zip(myLabels, q_xs, q_ys):
        matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                   textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                   bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    myAxes.set_title("Reciprocal lattice", y=1.05, fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-1, ymax=1, linewidth=0.5, color = 'b')
    myAxes.set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 22, rotation = 'vertical')
    myFigure.savefig("%s/intensitiesAtZero/measured_qx_qy_reciprocalLattice.png"%inputFolder, fontsize = 20)
    matplotlib.pyplot.close()
    
    # Plot of measured h, k points, with labels.
    myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
    myAxes = myFigure.add_subplot(1,1,1)    # Axes object  -  one row, one column, first plot    
    myAxes.scatter(hs_plot, ks_plot, color="red", marker=".")
    for label, x, y in zip(myLabels, hs_plot, ks_plot):
        matplotlib.pyplot.annotate(label, xy = (x, y), xytext = (-3, 3), size = 10, 
                                   textcoords = 'offset points', ha = 'right', va = 'bottom', 
                                   bbox = dict(boxstyle = 'round,pad=0.3', fc = 'yellow', alpha = 0.3, ec='none'))
    myAxes.set_title("Reciprocal lattice", y=1.05, fontsize = 20)
    myAxes.tick_params(axis='x', labelsize=14)
    myAxes.tick_params(axis='y', labelsize=14)
    matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
    matplotlib.pyplot.axvline(x=0, ymin=-1, ymax=1, linewidth=0.5, color = 'b')
    myAxes.set_xlabel("h", fontsize = 22, rotation = 'horizontal')
    myAxes.set_ylabel("k", fontsize = 22, rotation = 'vertical')
    myFigure.savefig("%s/intensitiesAtZero/measured_h_k_reciprocalLattice.png"%inputFolder, fontsize = 20)   
    matplotlib.pyplot.close()
    
if __name__ == "__main__":
    print "\n**** CALLING calculate_Patterson ****"
    calculate_PattersonFunction(sys.argv[1:])    