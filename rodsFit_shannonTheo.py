# -*- coding: utf-8 -*-
import matplotlib.pyplot
import getopt
import joblib 
import numpy
import scipy.optimize
import os
import sys

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
    

def shannonFit(myArguments):
    str_1 = '--resolutionLimit <resolutionLimit> --thickness <thickness>'
    str_2 = '--damping <damping> --folder <folder>'
    
    # DEFAULTS
    folder = './Output_runMergingVsModel'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["resolutionLimit=", 
                                               "thickness=", 
                                               "damping=", 
                                               "folder="])
    except getopt.GetoptError:
        print 'Usage: python rodsFit_shannonTheo.py %s'%(str_1, str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python rodsFit_shannonTheo.py %s'%(str_1, str_2)
            sys.exit()
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--thickness":
            d = float(value)
        elif option == "--damping":
            T = float(value)
        elif option == "--folder":
            folder = value           
    
    outputfolder = '%s/Shannon_sampling'%folder
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)  
    
    # DEFINE ROD INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
    
    shannon_model = []        
    print '%d Rods'%len(rodIndices)
    for rod in rodIndices:
        
        # EXTRACT ROD DATA
        rodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                 %(folder, rod[0], rod[1]))          
        experimental_q = rodObject.experimental_q
        experimental_I = rodObject.experimental_I
        qMin = rodObject.qMin
        qMax = rodObject.qMax
        
        # DEFINE SAMPLING POINTS
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
        N = (n-1)/2            # ALWAYS INTEGER
        params = [1 for i in range(0, n)]
        popt, pcov = scipy.optimize.curve_fit(lambda x, *params: 
                                              wrapper_fit_function(x, N, sampling, T, d, params), 
                                              experimental_q, experimental_I, p0=params)
        
        print "ROD ", rod
        print len(experimental_q), " EXPERIMENTAL POINTS"
        print "QROD EXTENSION (A-1): ", qMin, qMax
        print 'SAMPLING INTERVAL (A-1): ', sampling
        print 'SAMPLING AT: ', samplings                
        print 'N = ', N, ', n = ', n
        print 'FIT RESULTS: ', popt
        
        # GENERATE FITTING CURVE
        x = numpy.linspace(qMin, qMax, num=(qMax-qMin)/0.001, endpoint = True)
        #x = numpy.linspace(qMin, qMax, 10000)
        y = fit_function(x, popt, N, sampling, T, d)
        
        # MODEL
        for i in range(0, len(x)):
            spot_model = [int(rod[0]), int(rod[1]), x[i], y[i]] # h k qRod I
            shannon_model.append(spot_model)
                            
        # PLOT
        matplotlib.pyplot.figure()
        matplotlib.pyplot.scatter(experimental_q, experimental_I, color='c', alpha=0.1)
        matplotlib.pyplot.scatter(samplings, 
                                  [0 for i in range(0, len(samplings))], 
                                  color='m', 
                                  alpha=1)
        matplotlib.pyplot.plot(x, y, color='m')
        matplotlib.pyplot.savefig('%s/shannon_rod_%d_%d_nPts_%d.png'
                                   %(outputfolder, rod[0], rod[1], n), dpi=96*2)
        matplotlib.pyplot.close()
        
        # UPDATE BRAGG ROD WITH SHANNON MODEL COEFFICIENTS
        rodObject.setModelCoefficients(popt)
        if not os.path.exists("%s/braggRodObjects"%outputfolder):
            os.mkdir("%s/braggRodObjects"%outputfolder)
        joblib.dump(rodObject, '%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                %(outputfolder, rod[0], rod[1]))   
     
    # SAVE MODEL
    if not os.path.exists("%s/model"%outputfolder):
        os.mkdir("%s/model"%outputfolder)
    shannon_model = numpy.asarray(shannon_model, dtype=numpy.float32)
    joblib.dump(shannon_model, '%s/model/lattice_model.jbl'%outputfolder)
    
if __name__ == "__main__":
    print "\n**** CALLING ShannonFit ****"
    shannonFit(sys.argv[1:])    