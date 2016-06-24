# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 12:29:10 2016
@author: casadei_c
DRAW DISTRIBUTION OF N OF PEAKS (FROM peaks.txt FILE) PER IMAGE
"""
import matplotlib.pyplot
def nPeaksDistributionHistogramFunction(runNumber, nPeaksList):
    matplotlib.pyplot.close()
    maxNpeaks = max(nPeaksList)
    myBins = range(0, maxNpeaks+20, 10)
    matplotlib.pyplot.figure(facecolor = 'w')
    (n, bins, patches) = matplotlib.pyplot.hist(nPeaksList, bins = myBins)
    maxNimgsPlus = max(n) + 2
    maxNimgsPlus = int(maxNimgsPlus)
    yRange = range(0, maxNimgsPlus, maxNimgsPlus/5)
    matplotlib.pyplot.yticks(yRange)
    matplotlib.pyplot.title("Histogram")
    matplotlib.pyplot.ylabel("Number of images")
    matplotlib.pyplot.xlabel("Number of peaks")
    matplotlib.pyplot.savefig('./Output_r%s/ExtractExperimentalInfo/r%s_Histogram_Nimages_Npeaks.png'%(runNumber, runNumber), dpi = 300, facecolor = 'w')
    matplotlib.pyplot.close()
    print 'Saving histogram in ./Output_r%s/ExtractExperimentalInfo/r%s_Histogram_Nimages_Npeaks.png'%(runNumber, runNumber)