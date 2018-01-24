# -*- coding: utf-8 -*-
import numpy
import os
import pickle
import sys
import getopt
import matplotlib

import imageSums_utilities


            
def sumTilted_Function_spotWidth(myArguments):
    nCountsPerPhoton = 26
    
    str1 = '--halfWidth <halfWidth> --tiltAngle <tiltAngle>'
           
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["halfWidth=", 
                                                                 "tiltAngle="
                                                                 ])
    except getopt.GetoptError:
        print ('Error Usage: python imageSumming_tilted_spotWidth.py %s'
               %(str1))
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print ('Usage: python imageSumming_tilted_spotWidth.py %s'
                   %(str1))
            sys.exit()
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--tiltAngle":
            tiltAngle_deg = value
  
    # FOLDERS  
    inputFolder = './Output_imageSum_tilt_%s'%tiltAngle_deg  
    outputFolder = './Output_imageSum_gaussFit_tilt_%s'%tiltAngle_deg  
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
        
    sums_file = open('%s/sumMatrix_dictionary_list_tilt_%s.pkl'%(inputFolder, tiltAngle_deg), 'rb')
    dictionaryList = pickle.load(sums_file)
    sums_file.close()
    
    print len(dictionaryList)
    dictionaryList_gaussFit = []
    for dictionary in dictionaryList:
        h = dictionary['h']
        k = dictionary['k']
        qRod = dictionary['qRod']
        nTerms = dictionary['nTerms']
        sumMatrix = dictionary['sumMatrix']
        
        ### GAUSS FIT ###
        refined_sigma_x, refined_sigma_y, \
        refined_x0, refined_y0, \
        refined_amplitude, gauss_integral, \
        data, data_fitted = imageSums_utilities.do_gaussFit(sumMatrix)
        
        print h, k, qRod, gauss_integral/nCountsPerPhoton
        
        ### PLOT GAUSS FIT ###
        if not numpy.isnan(gauss_integral):
            myFigureObject = matplotlib.pyplot.figure()
            myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*halfWidth, 2*halfWidth), 
                                                         origin='lower', interpolation='nearest')
            try:
                matplotlib.pyplot.gca().contour(numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                numpy.linspace(0, 2*halfWidth-1, 2*halfWidth), 
                                                data_fitted.reshape(2*halfWidth, 2*halfWidth), 4, colors='w')
            except:
                print 'PROBLEM DRAWING CONTOURS'
                
            matplotlib.pyplot.title('Orbit: %d %d %.2f Gauss integral: %.1f counts\nSig_x %.2f Sig_y %.2f'
                                     %(h, k, qRod, gauss_integral, refined_sigma_x, refined_sigma_y))
            myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
            matplotlib.pyplot.savefig('%s/gauss_fit_%d_%d_%.2f_nTerms_%d.png'
                                       %(outputFolder, h, k, qRod, nTerms), dpi = 2*96)
            matplotlib.pyplot.close()  
        dictionary['refined_sigma_x'] = refined_sigma_x
        dictionary['refined_sigma_y'] = refined_sigma_y
        dictionary['refined_x0'] = refined_x0
        dictionary['refined_y0'] = refined_y0
        dictionary['refined_amplitude'] = refined_amplitude
        dictionary['gauss_integral'] = gauss_integral
        dictionaryList_gaussFit.append(dictionary)
        
    sums_file_gaussFit = open('%s/sumMatrix_dictionary_list_gaussFit_tilt_%s.pkl'%(outputFolder, tiltAngle_deg), 'wb')
    pickle.dump(dictionaryList_gaussFit, sums_file_gaussFit)
    sums_file_gaussFit.close()


if __name__ == "__main__":
    print "\n**** CALLING imageSumming_tilted_spotWidth ****"
    sumTilted_Function_spotWidth(sys.argv[1:])   