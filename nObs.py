# -*- coding: utf-8 -*-
import matplotlib.pyplot
import getopt
import joblib 
import numpy
import sys

import makeOrbits



def count_nObs(myArguments):
    
    str_1 = '--resolutionLimit <resolutionLimit> --cellSize <cellSize>'
    str_2 = '--bin_spacing <bin_spacing> --folder <folder>'
    
    # DEFAULTS
    bin_spacing = 0.01
    folder = './Output_runMergingVsModel'
        
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["resolutionLimit=", 
                                               "cellSize=", 
                                               "bin_spacing=", 
                                               "folder="])
    except getopt.GetoptError:
        print 'Usage: python nObs.py %s %s'%(str_1, str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python nObs.py %s %s'%(str_1, str_2)
            sys.exit()
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
        elif option == "--bin_spacing":
            bin_spacing = float(value)
        elif option == "--folder":
            folder = value  
    
    # RECIPROCAL BASIS VECTORS        
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],
                                          [0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I     
    
    # DEFINE ROD INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
    
    # LOGFILE 
    fOpen = open('%s/nObs.txt'%folder, 'w')
    fOpen.write('h k qRod q_3D d_3D N_obs')
    
    # CALCULATE N OBSERAVATIONS VS Q3D
    N_bins = []
    resolution_bins = []
     
    print '%d Rods'%len(rodIndices)
    for rod in rodIndices:
        
        # EXTRACT ROD DATA
        rodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                 %(folder, rod[0], rod[1]))          
        experimental_q = rodObject.experimental_q
        experimental_I = rodObject.experimental_I
        experimental_I_cleaned = [experimental_I[i] 
                                  for i in range(0, len(experimental_I)) 
                                  if not numpy.isnan(experimental_I[i])]
        experimental_q_cleaned = [experimental_q[i] 
                                  for i in range(0, len(experimental_I)) 
                                  if not numpy.isnan(experimental_I[i])]
        qMin = rodObject.qMin
        qMax = rodObject.qMax
        
        # DEFINE SAMPLING POINTS        
        bin_centers= []
        for n in range(0, 1000):
            bin_center = n * bin_spacing
            if bin_center <= qMax:
                bin_centers.append(bin_center)
            else:
                break        
        n = len(bin_centers)     # ONLY NON_NEGATIVE
        for i in range(1, n):
            bin_centers.append(-bin_centers[i])            
        bin_centers.sort()        
        n = len(bin_centers)     # ALWAYS ODD
        if not n%2 == 1:
            raise Exception('Even n.')
        print 'ROD ', rod
        print 'qRod range: [%.3f %.3f]'%(qMin, qMax)
        print 'Bin centers:'
        print bin_centers
        
        # COUNT N OBS PER BIN
        for bin_center in bin_centers:
            left  = bin_center - 1.5 * bin_spacing
            right = bin_center + 1.5 * bin_spacing
            Is_bin = [experimental_I_cleaned[i] 
                      for i in range(0, len(experimental_I_cleaned)) 
                      if (left < experimental_q_cleaned[i] < right)]
            N_bin = len(Is_bin)
            
            h = int(rod[0])
            k = int(rod[1])
            qRod = bin_center
            reciprocalVector = [h, k]*reciprocalCellRows
            q_x = reciprocalVector[0,0]         # A^(-1)
            q_y = reciprocalVector[0,1]         # A^(-1)
            q = numpy.sqrt(q_x**2 + q_y**2 + qRod**2)     # A^(-1)
            if q != 0:
                d_bin = 2* numpy.pi / q    # 
                
            N_bins.append(N_bin)
            resolution_bins.append(q)
            
            fOpen.write('%6d %6d %6.1f %6.1f %6.1f %8d \n'
                         %(h, k, qRod, q, d_bin, N_bin))
            
    fOpen.close()
    
    # PLOT BINNING
    plot_bin_width = 0.1
    plot_bins = numpy.arange(0, 2, plot_bin_width)
    plot_bin_qs = []
    plot_bin_Ns = []
    for i in range(0, len(plot_bins)-1):
        left = plot_bins[i]
        right = plot_bins[i+1]
        q_bin = (left+right)/2
        Ns_bin = [N_bins[j] for j in range(0, len(N_bins)) 
                            if (left < resolution_bins[j] < right)]
        if len(Ns_bin) > 0:
            N_bin = numpy.average(Ns_bin)
            plot_bin_qs.append(q_bin)
            plot_bin_Ns.append(N_bin)
    
    # PLOT N OBS IN ONE BIN VS Q3D
    matplotlib.pyplot.figure()
    title_str = 'N_obs in overlapping q_3D bins with spacing'
    matplotlib.pyplot.title(r'%s %.2f $\AA^{-1}$ and width %.2f $\AA^{-1}$'
                            %(title_str, bin_spacing, 3*bin_spacing), y=1.08)
    matplotlib.pyplot.scatter(resolution_bins, N_bins, color='c', alpha=0.1)
    matplotlib.pyplot.scatter(plot_bin_qs, 
                              plot_bin_Ns, 
                              s=120, facecolors='none', edgecolors='m')
    matplotlib.pyplot.gca().set_yscale('log')
    matplotlib.pyplot.gca().set_ylabel(r'$N_{\rm obs}$', fontsize = 17)
    matplotlib.pyplot.gca().set_xlabel(r'$q$ (3D) [$\AA^{-1}$]', fontsize = 17)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/nObs.png'%folder)

    
if __name__ == "__main__":
    print "\n**** CALLING nObs ****"
    count_nObs(sys.argv[1:])    