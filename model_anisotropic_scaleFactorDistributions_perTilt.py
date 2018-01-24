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
        print 'Usage: python model_anisotropic_scaleFactorDistributions_perTilt.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_anisotropic_scaleFactorDistributions_perTilt.py'
            sys.exit()
        
    tilts = ['5', '15', '20']
    runs_dictionary = {}
    runs_dictionary['5']  = ['0198', '0199']
    runs_dictionary['15'] = ['0195', '0196', '0197']
    runs_dictionary['20'] = ['0200', '0201']
            
    # CALCULATE NORMALIZED SCALES
    for tilt in tilts:
        print '******* TILT %s *******'%tilt
        nPos_Bperp = 0
        nNeg_Bperp = 0
        Ks = []
        B_paras = []
        B_perps = []
        runs = runs_dictionary['%s'%tilt]        
        for run in runs:           
            baseName = './Output_runMergingVsModel'
            inputFolder = '%s/transformAndScaleToModel_r%s'%(baseName, run)
            inputFolder_run = '%s/r%s-scalesVsModel_anisotropic_14A'%(inputFolder, 
                                                                  run)
            scales = joblib.load('%s/r%s-scalesVsModel_anisotropic.jbl'
                                 %(inputFolder_run, 
                                   run))   
              
            # GET SCALES
            print '******* RUN %s *******'%run          
            for lattice in range(0, len(scales)):
                factors = scales[lattice]
                scaleFactor  = factors[0]
                damping_para = factors[1]
                damping_perp = factors[2]
                
                if not numpy.isnan(scaleFactor):
                    Ks.append(scaleFactor)
                    B_paras.append(damping_para)
                    B_perps.append(damping_perp)
                    if damping_perp >= 0:
                        nPos_Bperp = nPos_Bperp + 1
                    else:
                        nNeg_Bperp = nNeg_Bperp + 1
            
        matplotlib.pyplot.figure(figsize=(30, 10))
        matplotlib.pyplot.title('B_perpendicular: %d positive, %d negative'%(nPos_Bperp,
                                                                             nNeg_Bperp),
                                                                             fontsize = 35)
        matplotlib.pyplot.hist(Ks, 
                               histtype='step', 
                               linewidth=3, 
                               label=r'$K_{\rm L}$')
        matplotlib.pyplot.hist(B_paras, 
                               histtype='step', 
                               linewidth=3, 
                               label=r'$B^{\parallel}_{\rm L}$')
        matplotlib.pyplot.hist(B_perps, 
                               histtype='step', 
                               linewidth=3, 
                               label=r'$B^{\perp}_{\rm L}$')
        matplotlib.pyplot.legend(frameon=False, fontsize=35)
        matplotlib.pyplot.gca().tick_params(axis='x', labelsize=25)
        matplotlib.pyplot.gca().tick_params(axis='y', labelsize=25)
        matplotlib.pyplot.savefig('%s/scale_distributions_%s_tilt_14A.png'%(baseName,
                                                                        tilt), 
                                                                        dpi=2*96)
        matplotlib.pyplot.close()
        print 'SAVEFIG'       
        
        
    
if __name__ == "__main__":
    print "\n**** CALLING model_anisotropic_scaleFactorDistributions ****"
    anisotropicScalesDistribution(sys.argv[1:])    