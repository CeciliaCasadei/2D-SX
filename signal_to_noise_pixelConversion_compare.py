# -*- coding: utf-8 -*-
import sys
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle
import warnings

###        
### IMAGE SUM - MODULE DISPLACE
###

def signalToNoise_compare(myArguments):
    warnings.filterwarnings("ignore")
 
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun="])
    except getopt.GetoptError:
        print 'python signalToNoise_pixelConversion_compare.py --selectedRun <selectedRun>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':            
            print 'python signalToNoise_pixelConversion_compare.py --selectedRun <selectedRun>'
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
       
    # FOLDERS
    outputFolder = './Output_r%s/Output_signal_to_noise_module_displace'%selectedRun
    
    colors = ['b', 'm', 'c']
    labels = ['', '_100_lattices', '_10_lattices']
    matplotlib.pyplot.figure()    
    
    StoN_compare = []
    
    index = 0
    for suffix in labels:
        
        # LOAD Qs
        fileToOpen = '%s/signalToNoise%s_qs.pkl'%(outputFolder, suffix)
        fRead = open(fileToOpen, 'rb')
        qs = pickle.load(fRead)                                               
        fRead.close()
        
        # LOAD S/N
        fileToOpen = '%s/signalToNoise%s_StoN.pkl'%(outputFolder, suffix)
        fRead = open(fileToOpen, 'rb')
        StoN = pickle.load(fRead)                                               
        fRead.close()
        
        # BINNING    
        q_plot = []
        signal_to_noise_plot = [] 
        
        bins = numpy.linspace(min(qs), max(qs), 25)
        for i in range(0, len(bins)-1):
            left_q  = bins[i]
            right_q = bins[i+1]
            signal_to_noise_bin = [StoN[i] for i in range(0, len(qs)) if left_q <= qs[i] <= right_q]
            q_mid_bin = (left_q+right_q)/2
            signal_to_noise_avg_bin = numpy.average(signal_to_noise_bin)
                
            q_plot.append(q_mid_bin)
            signal_to_noise_plot.append(signal_to_noise_avg_bin)

        matplotlib.pyplot.scatter(q_plot, signal_to_noise_plot, marker='o', s=6, color=colors[index])
                
        index = index + 1
        StoN_compare.append(signal_to_noise_plot)
        
    matplotlib.pyplot.gca().set_ylim([0.5, 1000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"$q$ ($\AA^{-1}$)", fontsize = 20, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"$S/N$", fontsize = 20, rotation = 'vertical')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_compare.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_compare.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    StoN_all = StoN_compare[0]
    StoN_100 = StoN_compare[1]
    StoN_10  = StoN_compare[2]
        
    ratios_all_to_10 = []
    ratios_100_to_10 = []
    ratios_all_to_100= []
    for idx in range(0, len(StoN_all)):
        ratio_all_to_10 = StoN_all[idx] / StoN_10[idx]
        ratio_100_to_10 = StoN_100[idx] / StoN_10[idx]
        ratio_all_to_100= StoN_all[idx] / StoN_100[idx]
        ratios_all_to_10.append(ratio_all_to_10)
        ratios_100_to_10.append(ratio_100_to_10)
        ratios_all_to_100.append(ratio_all_to_100)
            
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(q_plot, ratios_all_to_10)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_all_to_10.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_all_to_10.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(q_plot, ratios_100_to_10)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_100_to_10.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_100_to_10.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()
    matplotlib.pyplot.scatter(q_plot, ratios_all_to_100)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_all_to_100.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_ratios_all_to_100.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()


if __name__ == "__main__":
    print "\n**** CALLING signalToNoise_pixelConversion_compare ****"
    signalToNoise_compare(sys.argv[1:])