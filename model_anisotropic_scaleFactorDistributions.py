# -*- coding: utf-8 -*-
import getopt
import joblib
import numpy
import sys
import matplotlib.pyplot

def anisotropicScalesDistribution(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h")
    except getopt.GetoptError:
        print 'Usage: python model_applyScales_anisotropic.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyScales_anisotropic.py'
            sys.exit()
        
    runs = ['0199'] #, '0196', '0197', '0198', '0199', '0200', '0201']
            
    # CALCULATE NORMALIZED SCALES
            
    for run in runs:
        Ks = []
        B_paras = []
        B_perps = []
        outputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%run
        
        normalizedScales = joblib.load('%s/r%s-scalesVsModel_anisotropic/r%s-scalesVsModel_anisotropic.jbl'
                                       %(outputFolder, run, run))   
          
        # GET NORMALIZED SCALES
        print '******* RUN %s *******'%run          
        

        for lattice in range(0, len(normalizedScales)):
            factors = normalizedScales[lattice]
            scaleFactor  = factors[0]
            damping_para = factors[1]
            damping_perp = factors[2]
            
            if not numpy.isnan(scaleFactor):
                Ks.append(scaleFactor)
                B_paras.append(damping_para)
                B_perps.append(damping_perp)
        
        matplotlib.pyplot.figure(figsize=(30, 10))
        matplotlib.pyplot.hist(Ks, histtype='step', linewidth=3, label=r'$K_{\rm L}$')
        matplotlib.pyplot.hist(B_paras, histtype='step', linewidth=3, label=r'$B^{\parallel}_{\rm L}$')
        matplotlib.pyplot.hist(B_perps, histtype='step', linewidth=3, label=r'$B^{\perp}_{\rm L}$')
        matplotlib.pyplot.legend(frameon=False, fontsize=35)
        matplotlib.pyplot.gca().tick_params(axis='x', labelsize=25)
        matplotlib.pyplot.gca().tick_params(axis='y', labelsize=25)
        matplotlib.pyplot.savefig('%s/r%s_scale_distributions.png'%(outputFolder, run), dpi=2*96)
        matplotlib.pyplot.close()
        
        
        
    
if __name__ == "__main__":
    print "\n**** CALLING model_anisotropic_scaleFactorDistributions ****"
    anisotropicScalesDistribution(sys.argv[1:])    