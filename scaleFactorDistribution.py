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
        print 'Usage: python scaleFactorDistribution.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaleFactorDistribution.py --runNumber <runNumber>'
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
    matplotlib.pyplot.figure(figsize=(20, 20),facecolor = 'w')
    matplotlib.pyplot.rcParams['axes.linewidth'] = 3
    matplotlib.pyplot.tick_params(axis='both', which='major', length=18, width=3, labelsize=50, pad=20)
    matplotlib.pyplot.axvline(x=avgScale, color='r', linewidth=2)
    matplotlib.pyplot.axvline(x=avgScale-stdDev, color='r', linestyle=':', linewidth=3)
    matplotlib.pyplot.axvline(x=avgScale+stdDev, color='r', linestyle=':', linewidth=3)
    (n, bins, patches) = matplotlib.pyplot.hist(normalizedScales_clean, bins=myBins)
    matplotlib.pyplot.gca().text(0.97, 0.97, "$K_L=$ %.2f $\pm$ %.2f"%(avgScale, stdDev), horizontalalignment='right', verticalalignment='top', fontsize=48, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.xlabel(r'$K_L$', fontsize=60, labelpad=16)
    matplotlib.pyplot.ylabel(r"$n_{\rm lattices}$", fontsize=60)
    matplotlib.pyplot.savefig('./Output_r%s/transformAndScale/scaleFactorDistribution.png'%(runNumber), dpi = 300, facecolor = 'w')
    matplotlib.pyplot.savefig('./Output_r%s/transformAndScale/scaleFactorDistribution.pdf'%(runNumber), dpi = 300, facecolor = 'w')
    matplotlib.pyplot.close()
    
    print 'Scale factor distribution histogram plotted in: ./Output_r%s/transformAndScale/scalefactorDistribution.png'%(runNumber)
    

if __name__ == "__main__":
    print "\n**** CALLING scaleFactorDistribution ****"
    scaleFactorDistributionHistogramFunction(sys.argv[1:])    
    