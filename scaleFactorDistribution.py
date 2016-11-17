# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg') # Force matplotlib not to use any Xwindows backend.
import matplotlib.pyplot
import joblib
import numpy

def scaleFactorDistributionHistogramFunction(myArguments):
    runNumber = ''
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python cellSizeDistribution.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python cellSizeDistribution.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
    
    filename = './Output_r%s/transformAndScale/r%s-finalScales/r%s-finalScales-normalized.jbl'%(runNumber, runNumber, runNumber)
    normalizedScales = joblib.load(filename)    
    normalizedScales_clean = [normalizedScales[i] for i in range(0, len(normalizedScales)) if not numpy.isnan(normalizedScales[i])]
    print len(normalizedScales_clean), 'scaled, out of ', len(normalizedScales)
    normalizedScales_clean = numpy.ravel(normalizedScales_clean)
    
    avgScale = numpy.average(normalizedScales_clean)
    stdDev = numpy.std(normalizedScales_clean)
    print avgScale, stdDev
    
    matplotlib.pyplot.close()
    maxScale = max(normalizedScales_clean)
    minScale = min(normalizedScales_clean)
    myBins = numpy.linspace(minScale-0.1, maxScale+0.1, 15)
    matplotlib.pyplot.figure(facecolor = 'w')
    matplotlib.pyplot.axvline(x=avgScale, color='r')
    matplotlib.pyplot.axvline(x=avgScale-stdDev, color='r', linestyle=':')
    matplotlib.pyplot.axvline(x=avgScale+stdDev, color='r', linestyle=':')
    (n, bins, patches) = matplotlib.pyplot.hist(normalizedScales_clean, bins=myBins)
    matplotlib.pyplot.title("Scale factors: %.2f $\pm$ %.2f"%(avgScale, stdDev), y=1.04)
    matplotlib.pyplot.xlabel("Scale factor")
    matplotlib.pyplot.ylabel("Number of lattices")
    matplotlib.pyplot.savefig('./Output_r%s/transformAndScale/scaleFactorDistribution.png'%(runNumber), dpi = 300, facecolor = 'w')
    matplotlib.pyplot.close()
    
    print 'Scale factor distribution histogram plotted in: ./Output_r%s/transformAndScale/scalefactorDistribution.png'%(runNumber)
    

if __name__ == "__main__":
    print "\n**** CALLING scaleFactorDistribution ****"
    scaleFactorDistributionHistogramFunction(sys.argv[1:])    
    