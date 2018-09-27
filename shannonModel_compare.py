# -*- coding: utf-8 -*-
import matplotlib.pyplot
import getopt
import joblib 
import numpy
import os
import sys

import makeOrbits
import shannonSamplings

from bins import binLimits

def plot(myArguments):
    in_1 = '--resolutionLimit <resolutionLimit> --thickness <thickness>'
    in_2 = '--overSampling <overSampling>'
   
   
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["resolutionLimit=", 
                                               "thickness=",
                                               "overSampling="])
    except getopt.GetoptError:
        print 'Usage: python shannonModel_compare.py %s %s'%(in_1, in_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python shannonModel_compare.py %s %s'%(in_1, in_2)
            sys.exit()
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--thickness":
            d = float(value)
        elif option == "--overSampling":
            overSampling = value
         
    
    # DEFINE ROD INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
    
    f_Folder = './Output_runMergingVsModel'
    i_Folder = './Output_runMerging' 
    
    f_Model = joblib.load('%s/Shannon_sampling/model/lattice_model.jbl'%f_Folder)
    i_Model = joblib.load('%s/Shannon_sampling/model/lattice_model.jbl'%i_Folder)
   
    
    outputfolder = '%s/Shannon_sampling/Figures_OS_%s'%(f_Folder, overSampling)
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder) 
    
    
    print '%d Rods'%len(rodIndices)
    for rod in rodIndices:
        
        h = rod[0]
        k = rod[1]
        
        # EXTRACT ROD DATA
        rodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                 %(f_Folder, rod[0], rod[1]))          
        experimental_q = rodObject.experimental_q
        experimental_I = rodObject.experimental_I
        qMax = rodObject.qMax
        
        # DEFINE SAMPLING POINTS
        sampling = (2*numpy.pi) / (2*d)
        samplings = shannonSamplings.get_shannonSamplings(sampling, qMax)
        n = len(samplings)     # ALWAYS ODD
        if not n%2 == 1:
            raise Exception('Even n.')
        
        print "ROD ", rod
        print len(experimental_q), " EXPERIMENTAL POINTS"
        
        # EXTRACT INITIAL MODEL
        x_i = []
        y_i = []
        
        for i in i_Model:
            if i[0] == h and i[1] == k:
                x_i.append(i[2])
                y_i.append(i[3])
                
                
        # EXTRACT FINAL MODEL
        x_f = []
        y_f = []
        
        for f in f_Model:
            if f[0] == h and f[1] == k:
                x_f.append(f[2])
                y_f.append(f[3])
         
        # EXTRACT FW VALUES
        qs_FW    = []
        Is_FW    = []
        sigIs_FW = []
        
        nBins = len(binLimits) - 1
        for i in range(0, nBins):
            FW_data = joblib.load('%s/Shannon_sampling/French_Wilson_OS_%s/FW_uniques_bin_%d.jbl'
                                  %(f_Folder, overSampling, i))
            
            for spot in FW_data:
                
                h_spot  = int(spot[0])
                k_spot  = int(spot[1])
                                
                if h_spot == h and k_spot == k:
                    qRod    = float(spot[3])
                    I_FW    = float(spot[6])
                    sigI_FW = float(spot[7])
                
                    qs_FW.append(qRod)
                    Is_FW.append(I_FW)
                    sigIs_FW.append(sigI_FW)
            
            
                            
        # PLOT
        matplotlib.pyplot.figure()
        matplotlib.pyplot.scatter(experimental_q, 
                                  experimental_I, 
                                  color='xkcd:pale yellow', 
                                  alpha=0.1, 
                                  zorder=-1)
        matplotlib.pyplot.scatter(samplings, 
                                  [0 for i in range(0, len(samplings))], 
                                  color='k', 
                                  marker='x',
                                  alpha=1, 
                                  s=12)
        
        matplotlib.pyplot.plot(x_i, y_i, color='b')
        matplotlib.pyplot.plot(x_f, y_f, color='m', zorder=0)
        matplotlib.pyplot.scatter(qs_FW, Is_FW, color='k', s=1, zorder=1)
        matplotlib.pyplot.errorbar(qs_FW, 
                                   Is_FW, 
                                   yerr = sigIs_FW,
                                   linestyle="None",
                                   elinewidth=0.3,# width of error bar line
                                   ecolor='k',    # color of error bar
                                   capsize=3,     # cap length for error bar
                                   capthick=0.3   # cap thickness for error bar)
                                   )
        myAxis = matplotlib.pyplot.gca()
        matplotlib.pyplot.axhline(y=0, 
                                  xmin=-1, 
                                  xmax=1, 
                                  linewidth=0.5, 
                                  color = 'k')
        matplotlib.pyplot.axhline(y=10, 
                                  xmin=-1, 
                                  xmax=1, 
                                  linewidth=0.5, 
                                  color = 'b', 
                                  ls='dashed')
        myAxis.set_xlim([-0.60,+0.60])
        matplotlib.pyplot.xlabel(r'$q_{\rm rod}$ $[\AA^{-1}]$', fontsize=30)
        matplotlib.pyplot.ylabel('Intensity [photons]', fontsize=30)
        scale = 1.1*max(experimental_I)
        myAxis.set_ylim([-0.1*scale,1*scale])
        matplotlib.pyplot.gca().tick_params(axis='x', labelsize=18)
        matplotlib.pyplot.gca().tick_params(axis='y', labelsize=18)
        matplotlib.pyplot.tight_layout()
        matplotlib.pyplot.savefig('%s/y_Cmp_shannon_rod_%d_%d_nPts_%d.png'
                                   %(outputfolder, rod[0], rod[1], n), dpi=96*3)
        matplotlib.pyplot.close()
        
        
     
    
if __name__ == "__main__":
    print "\n**** CALLING plot ****"
    plot(sys.argv[1:])    