# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle



def signalToNoise_module_displace_plot():
    
    ###        
    ### IMAGE SUM - MODULE DISPLACE
    ###
    
    # PARAMETERS
    selectedRun = '0127'
 
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

    # BINNING    
    q_plot = []
    signal_to_noise_plot = [] 
    signal_to_noise_std_plot = [] 
    min_n_terms_plot = []
    min_n_terms_std_plot = []

    bins = numpy.linspace(min(q_list), max(q_list), 25)
    for i in range(0, len(bins)-1):
        left_q  = bins[i]
        right_q = bins[i+1]
        q_bin               = [q_list[i]               for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]    
        signal_to_noise_bin = [signal_to_noise_list[i] for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
        min_n_terms_bin     = [min_n_terms_list[i]     for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
        
        print len(q_bin)
        #q_avg_bin = numpy.average(q_bin)
        q_mid_bin = (left_q+right_q)/2
        signal_to_noise_avg_bin = numpy.average(signal_to_noise_bin)
        signal_to_noise_std_bin = numpy.std(signal_to_noise_bin)
        min_n_terms_avg_bin = numpy.average(min_n_terms_bin)
        min_n_terms_std_bin = numpy.std(min_n_terms_bin)

        q_plot.append(q_mid_bin)
        signal_to_noise_plot.append(signal_to_noise_avg_bin)
        signal_to_noise_std_plot.append(signal_to_noise_std_bin)
        min_n_terms_plot.append(min_n_terms_avg_bin)
        min_n_terms_std_plot.append(min_n_terms_std_bin)
        print q_mid_bin, signal_to_noise_avg_bin, signal_to_noise_std_bin, min_n_terms_avg_bin, min_n_terms_std_bin

    matplotlib.pyplot.figure()     
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=18, pad=4)    
    matplotlib.pyplot.scatter(q_list, signal_to_noise_list, s=2)
    matplotlib.pyplot.scatter(q_plot, signal_to_noise_plot, s=100, facecolors='c', edgecolors='none', alpha = 0.40)
    matplotlib.pyplot.gca().set_ylim([0.5, 1000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 27, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$S/N$", fontsize = 27, rotation = 'vertical')
    matplotlib.pyplot.gca().text(0.03, 0.97, "(a)", horizontalalignment='left', verticalalignment='top', fontsize=30, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.gca().tick_params(axis='both', which='major', labelsize=18, pad=4)    
    matplotlib.pyplot.scatter(q_list, min_n_terms_list, s=2)
    matplotlib.pyplot.scatter(q_plot, min_n_terms_plot, s=100, facecolors='c', edgecolors='none', alpha = 0.40)
    matplotlib.pyplot.gca().set_ylim([0.001, 10000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 27, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$N_{S/N=1}$", fontsize = 27, rotation = 'vertical')
    matplotlib.pyplot.gca().text(0.03, 0.97, "(b)", horizontalalignment='left', verticalalignment='top', fontsize=30, transform = matplotlib.pyplot.gca().transAxes) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/min_n_terms_module_displace.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()

if __name__ == "__main__":
    print "\n**** CALLING signalToNoise_module_displace ****"
    signalToNoise_module_displace_plot()