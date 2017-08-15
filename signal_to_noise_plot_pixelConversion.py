# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle
import scipy.optimize
import warnings

import imageSums_utilities



def signalToNoise_module_displace_plot(myArguments):
    
    str1 = "--selectedRun <selectedRun> --nBins <nBins>"
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "nBins="])
    except getopt.GetoptError:
        print 'Usage: python signal_to_noise_plot_pixelConversion.py %s'%str1
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python signal_to_noise_plot_pixelConversion.py %s'%str1
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
        elif option == "--nBins":
            nbins = int(value)
           
    ###        
    ### IMAGE SUM - MODULE DISPLACE
    ###
    
    warnings.filterwarnings("ignore")
 
    # FOLDERS
    outputFolder = './Output_r%s/Output_signal_to_noise_module_displace'%selectedRun
    
    # EXTRACT SIGNAL TO NOISE RESULTS    
    qs_file = open('%s/signalToNoise_qs.pkl'%outputFolder, 'rb')
    q_list = pickle.load(qs_file)
    qs_file.close()
    
    StoN_file = open('%s/signalToNoise_StoN.pkl'%outputFolder, 'rb')
    signal_to_noise_list = pickle.load(StoN_file)
    StoN_file.close()
    
    minN_file = open('%s/signalToNoise_minN.pkl'%outputFolder, 'rb')
    min_n_terms_list = pickle.load(minN_file)
    minN_file.close()

    print len(q_list), len(signal_to_noise_list), len(min_n_terms_list)#, len(signal_to_noise_k_list), len(min_n_terms_k_list)
    
    # BINNING    
    q_plot = []
    signal_to_noise_plot = [] 
    min_n_terms_plot = []

    bins = numpy.linspace(min(q_list), max(q_list), nbins)
    for i in range(0, len(bins)-1):
        left_q  = bins[i]
        right_q = bins[i+1]
        q_bin                     = [q_list[i]                     for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]    
        signal_to_noise_bin       = [signal_to_noise_list[i]       for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
        min_n_terms_bin           = [min_n_terms_list[i]           for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]

        print len(q_bin)
        q_mid_bin = (left_q+right_q)/2
        signal_to_noise_avg_bin = numpy.average(signal_to_noise_bin)
        min_n_terms_avg_bin = numpy.average(min_n_terms_bin)

        print q_mid_bin
        print min_n_terms_avg_bin
        print signal_to_noise_avg_bin
        print '\n'

        q_plot.append(q_mid_bin)
        signal_to_noise_plot.append(signal_to_noise_avg_bin)
        min_n_terms_plot.append(min_n_terms_avg_bin)

    matplotlib.pyplot.figure()     
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=6)    
    matplotlib.pyplot.scatter(q_list, signal_to_noise_list, s=2)
    matplotlib.pyplot.scatter(q_plot, signal_to_noise_plot, s=150, facecolors='c', edgecolors='none', alpha = 0.40)
    matplotlib.pyplot.gca().set_xlim([0, 1.7])
    matplotlib.pyplot.gca().set_ylim([0.5, 1000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 26, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$S/N$", fontsize = 26, rotation = 'vertical')
    matplotlib.pyplot.gca().text(0.03, 0.97, "(a)", horizontalalignment='left', verticalalignment='top', fontsize=30, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_%dbins.png'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_%dbins.pdf'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=16, pad=4)    
    matplotlib.pyplot.scatter(q_list, min_n_terms_list, s=2)
    matplotlib.pyplot.scatter(q_plot, min_n_terms_plot, s=150, facecolors='c', edgecolors='none', alpha = 0.40)
    matplotlib.pyplot.gca().set_xlim([0, 1.7])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 24, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$N_{S/N=1}$", fontsize = 24, rotation = 'vertical')
    matplotlib.pyplot.gca().text(0.03, 0.97, "(b)", horizontalalignment='left', verticalalignment='top', fontsize=30, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_%dbins.png'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_%dbins.pdf'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.close()

    q_plot_cleaned   = [q_plot[i]  for i in range(0, len(q_plot)) if ( q_plot[i] > 0.2 and not numpy.isnan(min_n_terms_plot[i]) )   ]
    min_n_cleaned    = [numpy.log(min_n_terms_plot[i])  for i in range(0, len(q_plot)) if ( q_plot[i] > 0.2 and not numpy.isnan(min_n_terms_plot[i]) )   ]
        
    popt, pcov = scipy.optimize.curve_fit(imageSums_utilities.line, q_plot_cleaned, min_n_cleaned)
    x = numpy.linspace(0.1, 1.6, 100)
    y_line = imageSums_utilities.line(x, *popt)
    I = numpy.exp(y_line)
        
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=22, pad=6)    
    matplotlib.pyplot.plot(x, I, 'm--')
    matplotlib.pyplot.scatter(q_list, min_n_terms_list, s=2)
    matplotlib.pyplot.scatter(q_plot, min_n_terms_plot, s=150, facecolors='c', edgecolors='none', alpha = 0.40)
    matplotlib.pyplot.gca().set_xlim([0, 1.7])
    matplotlib.pyplot.gca().set_ylim([0.001, 10000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 26, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$N_{S/N=1}$", fontsize = 26, rotation = 'vertical')
    matplotlib.pyplot.gca().text(0.03, 0.97, "(b)", horizontalalignment='left', verticalalignment='top', fontsize=30, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_%dbins_fit.png'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_%dbins_fit.pdf'%(outputFolder, nbins), dpi=4*96)
    matplotlib.pyplot.close()



if __name__ == "__main__":
    print "\n**** CALLING signal_to_noise_plot_pixelConversion ****"
    signalToNoise_module_displace_plot(sys.argv[1:])