# -*- coding: utf-8 -*
import numpy
import getopt
import sys
from scipy import interpolate
import matplotlib
from matplotlib import cm
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


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
            
    structureFactorsFile = open('./theo_ACF/1fbb_SF.txt', 'r')
           
    IsAtZero = []
    
    for structureFactor in structureFactorsFile:
        splitLine = structureFactor.split()
        h = int(splitLine[0])
        k = int(splitLine[1])
        F = float(splitLine[3])
        I = (F**2)
        print '%d %d %f %f'%(h, k, F, I)
        IatZero = []
        IatZero.append(h)
        IatZero.append(k)
        IatZero.append(I)
        IsAtZero.append(IatZero)
    IsAtZero = numpy.asarray(IsAtZero)
   
    # Matrix: h k I(qRod = 0)                    
    totalI = 0
    for IatZero in IsAtZero:
        intensity = IatZero[2]
        totalI = totalI + intensity
    averageI = totalI/IsAtZero.shape[0]
    
    # Normalized matrix: h k I(qRod = 0)
    normalizedIsAtZero = []
    for IatZero in IsAtZero:
        intensity = IatZero[2]
        normalizedIntensity = intensity/averageI
        normalizedIatZero = [IatZero[0], IatZero[1], normalizedIntensity]
        normalizedIsAtZero.append(normalizedIatZero)
    normalizedIsAtZero = numpy.asarray(normalizedIsAtZero)
    
    print IsAtZero 
    print normalizedIsAtZero
    
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
            for h in h_s:      # from -8 to +8
                for k in k_s:  # from -8 to +8
                    
                    # Match h, k with corresponding Bragg rod.
                    nMatches = 0
                    for braggRod in normalizedIsAtZero:
                        hRod = braggRod[0]
                        kRod = braggRod[1]
                        if (   (h == hRod and k == kRod) 
                            or (h == -hRod-kRod and k == hRod) 
                            or (h == kRod and k == -hRod-kRod) 
                            or (h == -hRod and k == -kRod) 
                            or (h == hRod+kRod and k == -hRod) 
                            or (h == -kRod and k == hRod+kRod)   ):
                                
                            pattersonSum = pattersonSum + braggRod[2] * numpy.cos(2*numpy.pi*(h*u + k*v))
                            nMatches = nMatches + 1
                            
                            ### PLOT ###
                            hs_plot.append(h)                   # Measured h
                            ks_plot.append(k)                   # Measured k
                            reciprocalVector = [h, k]*reciprocalCellRows
                            q_x = reciprocalVector[0,0]         # A^(-1)
                            q_y = reciprocalVector[0,1]         # A^(-1)                               
                            q_xs.append(q_x)                    # Measured qx
                            q_ys.append(q_y)                    # Measured qy
                            Is_plot.append(braggRod[2])         # Corresponding measured I(qRod = 0), imposed P6
                            ###
                            
                            if nMatches > 1:
                                print 'nMatches > 1'
                                
            xCart = u * a + v * a * numpy.cos(2*numpy.pi/3)           
            yCart = v * a * numpy.sin(2*numpy.pi/3)
            xCart_s.append(xCart)
            yCart_s.append(yCart)
            patterson.append(pattersonSum)
            
            
    ############FIGURES###################
    print min(q_xs)
    print max(q_xs)
    print min(q_ys)
    print max(q_ys)
    
    # 2D plots of measured h, k; measured qx, qy and x, y point where patterson was calculated
    myFigure = matplotlib.pyplot.figure(figsize=(15, 15), facecolor = 'w')     # Figure object
    
    myAxes_1 = myFigure.add_subplot(3,1,1)    # Axes object  -  one row, one column, first plot    
    myAxes_1.scatter(hs_plot, ks_plot, color="red", marker=".")
    myAxes_1.tick_params(axis='x', labelsize=14)
    myAxes_1.tick_params(axis='y', labelsize=14)
    myAxes_1.set_xlabel("h", fontsize = 18, rotation = 'horizontal')
    myAxes_1.set_ylabel("k", fontsize = 18, rotation = 'vertical')
    myAxes_1.axis('equal')
    
    myAxes_2 = myFigure.add_subplot(3,1,2)    # Axes object  -  one row, one column, first plot    
    myAxes_2.scatter(q_xs, q_ys, color="red", marker=".")
    myAxes_2.tick_params(axis='x', labelsize=14)
    myAxes_2.tick_params(axis='y', labelsize=14)
    myAxes_2.set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 18, rotation = 'horizontal')
    myAxes_2.set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 18, rotation = 'vertical')
    myAxes_2.axis('equal')
    
    myAxes_3 = myFigure.add_subplot(3,1,3)    # Axes object  -  one row, one column, first plot    
    myAxes_3.scatter(xCart_s, yCart_s, color="red", marker=".")
    myAxes_3.tick_params(axis='x', labelsize=14)
    myAxes_3.tick_params(axis='y', labelsize=14)
    myAxes_3.set_xlabel("x", fontsize = 18, rotation = 'horizontal')
    myAxes_3.set_ylabel("y", fontsize = 18, rotation = 'vertical')
    myAxes_3.axis('equal')
        
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/1fbb_real_reciprocalSpace.png'%inputFolder)
    matplotlib.pyplot.close()
    
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
    myFigure.savefig("%s/intensitiesAtZero/1fbb_measured_qx_qy_reciprocalLattice.png"%inputFolder, fontsize = 20)
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
    myFigure.savefig("%s/intensitiesAtZero/1fbb_measured_h_k_reciprocalLattice.png"%inputFolder, fontsize = 20)   
    matplotlib.pyplot.close()
    
    # 3D plot of I(qx, qy) at qRod = 0        
    myFigure = matplotlib.pyplot.figure()
    myAxes = myFigure.gca(projection='3d')
    mySurface = myAxes.plot_trisurf(q_xs, q_ys, Is_plot, cmap=cm.seismic, linewidth=0.2)
    myAxes.set_xticklabels([-1, '', 0, '', 1])
    myAxes.set_yticklabels([-1, '', 0, '', 1])
    myAxes.set_zticklabels([])
    myFigure.colorbar(mySurface, shrink=0.5, aspect=5)
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    myFigure.tight_layout()
    matplotlib.pyplot.show()
    matplotlib.pyplot.close()
    
    # Interpolated plot of I(qx, qy) at qRod = 0  
    qx_s_interpolate = numpy.linspace(-1.0, +1.0, 100)
    qy_s_interpolate = numpy.linspace(-1.0, +1.0, 100)    
    (q_xx, q_yy) = numpy.meshgrid(qx_s_interpolate, qy_s_interpolate) # Interpolation grid: 100x100, 100x100
    
    inputCoos = numpy.zeros(shape=(len(q_xs), 2))   # 210x2
    for i in range(0, inputCoos.shape[0]):
        inputCoos[i, 0] = q_xs[i]
        inputCoos[i, 1] = q_ys[i]
        
    f = interpolate.griddata(inputCoos, Is_plot, (q_xx, q_yy), method='linear', fill_value=numpy.nan )  # 100x100
    
    myAxesImageObject = matplotlib.pyplot.imshow(f, origin='lower', interpolation='nearest', extent=[-1.0, 1.0, -1.0, 1.0], vmin = 0.3, vmax = 1.2)
    matplotlib.pyplot.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.gca().set_xlabel("q$_x$ (A$^{-1}$)", fontsize = 15, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel("q$_y$ (A$^{-1}$)", fontsize = 15, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/1fbb_normalized_interpolated_IsAtZero.png'%inputFolder)
    matplotlib.pyplot.close()         
                      
    # 3D plot of patterson function vs real space        
    fig = matplotlib.pyplot.figure()
    ax = fig.gca(projection='3d')
    mySurface = ax.plot_trisurf(xCart_s, yCart_s, patterson, cmap=cm.seismic, linewidth=0.2)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    fig.colorbar(mySurface, shrink=0.5, aspect=5)
    matplotlib.pyplot.axis('off')
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    matplotlib.pyplot.show()
    
    # Interpolated plot of patterson function vs real space
    x_s_interpolate = numpy.linspace(-70.0, +70.0, 140)
    y_s_interpolate = numpy.linspace(-70.0, +70.0, 140)    
    (xx, yy) = numpy.meshgrid(x_s_interpolate, y_s_interpolate) 
    
    inputCoos = numpy.zeros(shape=(len(xCart_s), 2))  
    for i in range(0, inputCoos.shape[0]):
        inputCoos[i, 0] = xCart_s[i]
        inputCoos[i, 1] = yCart_s[i]
        
    f = interpolate.griddata(inputCoos, patterson, (xx, yy), method='linear', fill_value=numpy.nan )  
    
    myAxesImageObject = matplotlib.pyplot.imshow(f, origin='lower', interpolation='nearest', extent=[-70.0, 70.0, -70.0, 70.0])
    matplotlib.pyplot.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.gca().set_xlabel("x (A)", fontsize = 15, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel("y (A)", fontsize = 15, rotation = 'vertical')
    matplotlib.pyplot.xlim([-40, 65])
    matplotlib.pyplot.ylim([-2, 55])
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/1fbb_normalized_interpolated_patterson.png'%inputFolder)
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
    mySurface = ax.plot_trisurf(xCart_s_extended, yCart_s_extended, patterson_extended, cmap=cm.seismic, linewidth=0.2)
    ax.set_zlim(-4500, 4500)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    fig.colorbar(mySurface, shrink=0.5, aspect=5)
    matplotlib.pyplot.axis('off')
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    matplotlib.pyplot.show()  
    


    # Interpolated plot of translation-extended patterson function vs real space
    x_s_interpolate = numpy.linspace(-70.0, +130.0, 200)
    y_s_interpolate = numpy.linspace(-70.0, +130.0, 200)    
    (xx, yy) = numpy.meshgrid(x_s_interpolate, y_s_interpolate) 
    
    inputCoos = numpy.zeros(shape=(len(xCart_s_extended), 2))  
    for i in range(0, inputCoos.shape[0]):
        inputCoos[i, 0] = xCart_s_extended[i]
        inputCoos[i, 1] = yCart_s_extended[i]
        
    f = interpolate.griddata(inputCoos, patterson_extended, (xx, yy), method='linear', fill_value=numpy.nan ) 
    
    myAxesImageObject = matplotlib.pyplot.imshow(f, origin='lower', interpolation='nearest', extent=[-70.0, 130.0, -70.0, 130.0])
    matplotlib.pyplot.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    matplotlib.pyplot.gca().set_xlabel("x (A)", fontsize = 15, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel("y (A)", fontsize = 15, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/1fbb_normalized_interpolated_patterson_translated.png'%inputFolder)
    matplotlib.pyplot.close()         

    
    # Plot unique sector I(qx, qy) with qRod = 0
    uniqueSector = []
    for normalizedIatZero in normalizedIsAtZero:
        h = normalizedIatZero[0]
        k = normalizedIatZero[1]
        I = normalizedIatZero[2]
        reciprocalVector = [h, k]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1) 
        uniqueSectorRow = [q_x, q_y, I]
        uniqueSector.append(uniqueSectorRow)
    uniqueSector = numpy.asarray(uniqueSector)
    
    fig = matplotlib.pyplot.figure()
    ax = fig.gca(projection='3d')
    mySurface = ax.plot_trisurf(uniqueSector[:, 0], uniqueSector[:, 1], uniqueSector[:, 2], cmap=cm.seismic, linewidth=0.2)
    ax.set_zticklabels([])
    fig.colorbar(mySurface, shrink=0.5, aspect=5)
    matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
    matplotlib.pyplot.show()
    
    # Interpolated plot of unique sector of I(qx, qy) at qRod = 0  
    qx_s_interpolate = numpy.linspace(-1.0, +1.0, 100)
    qy_s_interpolate = numpy.linspace(-1.0, +1.0, 100)
    (q_xx, q_yy) = numpy.meshgrid(qx_s_interpolate, qy_s_interpolate)
    
    inputCoos = numpy.zeros(shape=(uniqueSector.shape[0], 2))
    for i in range(0, inputCoos.shape[0]):
        inputCoos[i, 0] = uniqueSector[i, 0]
        inputCoos[i, 1] = uniqueSector[i, 1]

    f = interpolate.griddata(inputCoos, uniqueSector[:, 2], (q_xx, q_yy), method='linear', fill_value=numpy.nan )
    
    myAxesImageObject = matplotlib.pyplot.imshow(f, origin='lower', interpolation='nearest', extent=[-1.0, 1.0, -1.0, 1.0])
    matplotlib.pyplot.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
    
    matplotlib.pyplot.savefig('%s/intensitiesAtZero/1fbb_normalized_interpolated_unique_IsAtZero.png'%inputFolder)
    matplotlib.pyplot.close()
    
    
if __name__ == "__main__":
    print "\n**** CALLING calculate_Patterson ****"
    calculate_PattersonFunction(sys.argv[1:])    