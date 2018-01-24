# -*- coding: utf-8 -*-
import os
import numpy
import pickle
import sys
import getopt
import matplotlib.pyplot
import scipy.optimize

import imageSums_utilities
import makeOrbits

def fit_function(x, a, N, delta, T, d):
    y = 0
    for i in range(-N, N+1):
        j = i+N
        y = (y + 
             a[j]  * (numpy.sin(d*(x-i*delta)) / (d*(x-i*delta)) ) 
                   * numpy.exp(-0.5*(1/(4*numpy.pi**2))*T*(x-i*delta)**2))
    return y
    
def wrapper_fit_function(x, N, delta, T, d, *args):
    n = 2*N + 1
    a = list(args[0][:n])
    return fit_function(x, a, N, delta, T, d)
            
def sumTilted_Function_finalIntegration(myArguments):
    
    str1 = '--ellipse_multiplicative_factor <ellipse_multiplicative_factor>'
    str2 = '--tiltAngle <tiltAngle> --resolutionLimit <resolutionLimit>'
    str3 = '--nominalCell <nominalCell> --thickness <thickness> --damping <damping>'
           
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["ellipse_multiplicative_factor=", 
                                                                 "tiltAngle=",
                                                                 "resolutionLimit=",
                                                                 "nominalCell=",
                                                                 "thickness=", 
                                                                 "damping="])
                                                                 
    except getopt.GetoptError:
        print ('Error Usage: python imageSumming_tilted_finalIntegration.py %s %s %s'
               %(str1, str2, str3))
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print ('Usage: python imageSumming_tilted_finalIntegration.py %s %s %s'
                   %(str1, str2, str3))
            sys.exit()
        elif option == "--ellipse_multiplicative_factor":
            ellipse_multiplicative_factor = float(value)
        elif option == "--tiltAngle":
            tiltAngle_deg = value
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--nominalCell":
            nominalCell = float(value)
        elif option == "--thickness":
            d = float(value)
        elif option == "--damping":
            T = float(value)
    
    # FOLDERS  
    inputFolder = './Output_imageSum_gaussFit_tilt_%s'%tiltAngle_deg  
    outputFolder = './Output_imageSum_gaussFit_tilt_%s/BraggLines'%tiltAngle_deg 
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)

    # RECIPROCAL CELL
    directCell = nominalCell * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],
                                             [0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2* numpy.pi * directCell.I  
    
    # LOAD PEAK SUMS    
    sums_file = open('%s/sumMatrix_dictionary_list_gaussFit_tilt_%s.pkl'
                     %(inputFolder, tiltAngle_deg), 'rb')
    dictionaryList = pickle.load(sums_file)
    sums_file.close()    
    print '%d peak sums.'%len(dictionaryList)
    
    fRead = open('%s/sigmaXCurveParameters.pkl'%inputFolder, 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()

    fRead = open('%s/sigmaYCurveParameters.pkl'%inputFolder, 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close() 
    
    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)

    
    # EXTRACT PEAK SUM
    for indices in rodIndices:
        print indices
        h_label = indices[0]
        k_label = indices[1]
        qRods = []
        Is = []
        for dictionary in dictionaryList:
            h = dictionary['h']
            k = dictionary['k']
            if h == h_label and k == k_label:
                
                qRod = dictionary['qRod']
                if qRod == 0:
                    qRod = qRod + 0.000000001
                bgSubtracted_total_sum = dictionary['sumMatrix'] 
                
                reciprocalVector = [h, k]*reciprocalCellRows
                q_x = reciprocalVector[0,0]         # A^(-1)
                q_y = reciprocalVector[0,1]         # A^(-1)
                q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1) 
                
                
                
                sigmaX = imageSums_utilities.poly_2_2D((q_2D, qRod), sigmaXCurveParameters[0], 
                                                                     sigmaXCurveParameters[1], 
                                                                     sigmaXCurveParameters[2])  
                                                                  
                sigmaY = imageSums_utilities.poly_2_2D((q_2D, qRod), sigmaYCurveParameters[0], 
                                                                     sigmaYCurveParameters[1], 
                                                                     sigmaYCurveParameters[2])  
                                                                     
                ### SUM INTEGRATION ON VARIABLE ELLIPSE ###
                integratedIntensity_ellipse = imageSums_utilities.integrate_ellipse(bgSubtracted_total_sum, 
                                                                                    sigmaX, 
                                                                                    sigmaY, 
                                                                                    ellipse_multiplicative_factor)
                qRods.append(qRod)
                Is.append(integratedIntensity_ellipse)
        
        
                
        # DEFINE SAMPLING POINTS
        qMax = max([max(qRods), -min(qRods)])
        sampling = (2*numpy.pi) / (2*d)
        samplings = []
        for n in range(0, 1000):
            q_rod_sampling = n * sampling
            if q_rod_sampling <= qMax:
                samplings.append(q_rod_sampling)
            else:
                samplings.append(q_rod_sampling)  # ADD ONE MORE TO BETTER FIT SHORT RODS
                break        
        n = len(samplings)     # ONLY NON_NEGATIVE
        for i in range(1, n):
            samplings.append(-samplings[i])            
        samplings.sort()        
        n = len(samplings)     # ALWAYS ODD
        if not n%2 == 1:
            raise Exception('Even n.')
        
        # FIT
        qRods = numpy.asarray(qRods)
        Is = numpy.asarray(Is)
        N = (n-1)/2            # ALWAYS INTEGER
        params = [1 for i in range(0, n)]
        popt, pcov = scipy.optimize.curve_fit(lambda x, *params: 
                                              wrapper_fit_function(x, N, sampling, T, d, params), 
                                              qRods, Is, p0=params)
        
        print "ROD ", h_label, k_label
        print len(qRods), " EXPERIMENTAL POINTS"
        print "QROD EXTENSION (A-1): ", -qMax, qMax
        print 'SAMPLING INTERVAL (A-1): ', sampling
        print 'SAMPLING AT: ', samplings                
        print 'N = ', N, ', n = ', n
        print 'FIT RESULTS: ', popt
        
        # GENERATE FITTING CURVE
        x = numpy.linspace(min(qRods), max(qRods), num=(2*qMax)/0.001, endpoint = True)
        y = fit_function(x, popt, N, sampling, T, d)
        
        # PLOT BRAGG LINES
        matplotlib.pyplot.figure()  
        matplotlib.pyplot.scatter(qRods, Is, edgecolors='none')
        matplotlib.pyplot.plot(x, y, color='m')
        matplotlib.pyplot.scatter(samplings, 
                                  [0 for i in range(0, len(samplings))], 
                                  color='m', 
                                  alpha=1)
        myAxis = matplotlib.pyplot.gca()
        matplotlib.pyplot.axhline(y=0, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
        matplotlib.pyplot.axhline(y=10, xmin=-1, xmax=1, linewidth=0.5, color = 'b')
        myAxis.set_xlim([-0.60,+0.60])
        scale = 1.1*max(Is)
        myAxis.set_ylim([-0.1*scale,1*scale])
        myAxis.set_xlabel("q$_z$ (A$^{-1}$)", fontsize = 12, rotation = 'horizontal')
        matplotlib.pyplot.savefig('%s/BraggLine_%d_%d'%(outputFolder, h_label, k_label))
        matplotlib.pyplot.close()
    
                   
   
        

if __name__ == "__main__":
    print "\n**** CALLING imageSumming_tilted_finalIntegration ****"
    sumTilted_Function_finalIntegration(sys.argv[1:])   