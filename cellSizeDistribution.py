# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg') # Force matplotlib not to use any Xwindows backend.
import matplotlib.pyplot
import joblib
import numpy

def cellSizeDistributionHistogramFunction(myArguments):
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
    
    filename = './Output_r%s/UnassembledImageProcessing/refinedLatticeSizes_r%s.jbl'%(runNumber, runNumber)
    cellSizes = joblib.load(filename)
    print len(cellSizes)
    
    avgCell = numpy.average(cellSizes)
    stdDev = numpy.std(cellSizes)
    
    matplotlib.pyplot.close()
    maxSize = max(cellSizes)
    minSize = min(cellSizes)
    myBins = numpy.linspace(minSize-0.05, maxSize+0.1, 20)
    matplotlib.pyplot.figure(figsize=(20, 20), facecolor = 'w')
    matplotlib.pyplot.tick_params(axis='both', which='major', length=18, width=3, labelsize=50, pad=20)
    (n, bins, patches) = matplotlib.pyplot.hist(cellSizes, bins=myBins)
    matplotlib.pyplot.gca().text(0.99, 0.98, "$a=$ (%.2f $\pm$ %.2f) $\AA$"%(avgCell, stdDev), horizontalalignment='right', verticalalignment='top', fontsize=48, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.xlabel(r'$a$ ($\AA$)', fontsize=60, labelpad=23)
    matplotlib.pyplot.ylabel(r"$n_{\rm lattices}$", fontsize=60)
    matplotlib.pyplot.savefig('./Output_r%s/cellSizeDistribution.png'%(runNumber), dpi = 300, facecolor = 'w')
    matplotlib.pyplot.close()
    
    print 'Cell size distribution histogram plotted in: ./Output_r%s/cellSizeDistribution.png'%(runNumber)
    

if __name__ == "__main__":
    print "\n**** CALLING cellSizeDistribution ****"
    cellSizeDistributionHistogramFunction(sys.argv[1:])    
    