# -*- coding: utf-8 -*-
import os
import sys
import getopt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy
import pickle
import joblib
import h5py



import imageSums_utilities


def signalToNoise_module_displace(myArguments):
    
    ###        
    ### IMAGE SUM - MODULE DISPLACE
    ###
    
    str1 = "--selectedRun <selectedRun> --nCountsPerPhoton <nCountsPerPhoton> --ellipse_multiplicative_factor <ellipse_multiplicative_factor> --precisionFactor <precisionFactor> --halfWidth <halfWidth> --label <label>"
    # READ INPUTS
    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun=", 
                                                                 "nCountsPerPhoton=", 
                                                                 "ellipse_multiplicative_factor=",
                                                                 "precisionFactor=",
                                                                 "halfWidth=",
                                                                 "label="]) # "!", "_10_lattices", "_100_lattices"
    except getopt.GetoptError:
        print 'Usage: python imageSums_displaced_modules_finalIntegration_pixelConversion.py %s'%str1
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python imageSums_displaced_modules_finalIntegration_pixelConversion.py %s'%str1
            sys.exit()
        elif option == "--selectedRun":
            selectedRun = value.zfill(4)
        elif option == "--nCountsPerPhoton":
            nCountsPerPhoton = int(value)
        elif option == "--ellipse_multiplicative_factor":
            ellipse_multiplicative_factor = float(value)
        elif option == "--precisionFactor":
            precision_factor = int(value)
        elif option == "--halfWidth":
            halfWidth = int(value)
        elif option == "--label":
            folder_label = value
    
    print "label: ", folder_label 
    folder_label = folder_label.strip("!")
    print "label: ", folder_label  
    
    k_squared = (10**precision_factor)**2 #100
  
    # FOLDERS
    outputFolder = './Output_r%s/Output_signal_to_noise_module_displace'%selectedRun
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)       
    
    # EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
    fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaXCurveParameters.pkl'%selectedRun, 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()
    
    fRead = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/sigmaYCurveParameters.pkl'%selectedRun, 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close()   
    
    # LOAD ORBITS
    fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements%s/orbits.pkl'%(selectedRun, folder_label)
    fRead = open(fileToOpen, 'rb')
    orbitsDictionary = pickle.load(fRead)                                               
    fRead.close()
    
    # PREPARE LATTICES TO IMAGES MATCHING
    imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%selectedRun    
    lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                           %(selectedRun, selectedRun, selectedRun)) # 629 lattices  
    lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(selectedRun, selectedRun, selectedRun), 'r')
    lattices_names = list(lattices_names)    
    images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(selectedRun, selectedRun))
    images_names = list(images_names)
    
    # LOG
    flog = open('%s/log%s.txt'%(outputFolder, folder_label), 'w')
    flog.write('#h, k, q, n_terms, n_ellipse, n_ring, I_img_ellipse_ph, avg_var_ring_pix_img, n_bgPixels\n')
    
    # RESULTS    
    q_list = []
    signal_to_noise_list = []
    min_n_terms_list = []
    
    # LOOP ON ORBITS
    for key, orbit in orbitsDictionary.items():        
        label = orbit.label        
        
        # h, k
        h_label = label[0]
        k_label = label[1]
        
        indices_1 = orbit.orbitIndices[0]
        indices_2 = orbit.orbitIndices[1]
        indices_3 = orbit.orbitIndices[2]
        h_1 = indices_1[0]
        k_1 = indices_1[1]
        h_2 = indices_2[0]
        k_2 = indices_2[1]
        h_3 = indices_3[0]
        k_3 = indices_3[1]
        
        # q
        q = 2 * numpy.pi / orbit.resolution
        
        # n_terms
        n_terms = orbit.nTerms
       
        ###
        # EXTRACT BG SUBTRACTED IMAGE SUM 
        # CONVERT TO PHOTON NUMBER 
        # KEEP NORMALIZATION ON N_TERMS
        ###
        bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtracted_total_sum # (300, 300) float64
        bgSubtractedOrbitSumMatrix_rotated = bgSubtractedOrbitSumMatrix_rotated / nCountsPerPhoton   
        
        # n_bgPixels used in bg plane calculation
        n_bgPixels = orbit.n_bgPixels / k_squared         # n of original pixels
        
        ###
        # BUILD MASKS 300x300
        ###
        sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
        integrationMask, ringMask = imageSums_utilities.buildMasks(bgSubtractedOrbitSumMatrix_rotated, 
                                                                   precision_factor, 
                                                                   ellipse_multiplicative_factor, 
                                                                   sigma_x, 
                                                                   sigma_y)
        
        # Total intensity on the ellipse, in the sum, expressed in photons, per image
        I_img_ellipse_ph = numpy.multiply(integrationMask, bgSubtractedOrbitSumMatrix_rotated).sum() / k_squared  # real I
        n_ellipse = (integrationMask.sum()) / k_squared   # n of original pixels
        
        ###
        # CALCULATE THE AVERAGE VARIANCE OF RING INTENSITY FROM SINGLE IMAGES avg_var_ring_pix_img
        ###
        var_ring_pix_img_list = []
        n_max = 100
        
        # Loop on lattices
        for index in range(0, len(lattices_names)):
            
            if len(var_ring_pix_img_list) > n_max:
                break
            
            lattice_name = lattices_names[index]
            imageNumber = int(lattice_name[80:84])
            latticeMatrix = lattices[index]                    # h_t k_t q_rod I_scaled flag i_unassembled j_unassembled scale
            image_name = images_names[imageNumber-1].split()[1]
            
            if latticeMatrix[0, 4] == 0:                       # Check transformation and scaling flag
                continue
           
            # LOAD UNASSEMBLED IMAGE
            unassembledDataFile = h5py.File('%s/%s'%(imagesDirectoryName, image_name), 'r')
            unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
            unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####  (1480, 1552)
            nRows = unassembledData.shape[0]
            nColumns = unassembledData.shape[1]
            
            # IN THE CURRENT LATTICE, LOOP ON SPOTS BELONGING TO THE CURRENT ORBIT  
            for spot in latticeMatrix:
                h = int(spot[0])
                k = int(spot[1])
                if (h == h_1 and k == k_1) or (h == h_2 and k == k_2) or (h == h_3 and k == k_3):
                    i = int(spot[5])
                    j = int(spot[6])                   
                    
                    left_edge   = j - halfWidth
                    right_edge  = j + halfWidth
                    bottom_edge = i - halfWidth
                    top_edge    = i + halfWidth
                    if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                        continue
                    
                    ### EXTRACT SECTOR ###
                    spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]  #(50, 50) float32
                    
                    # CONVERT SINGLE IMGS IN N OF PH
                    spotMatrix = spotMatrix / nCountsPerPhoton
                    
                    ### SCALE ###
                    spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                    spotMatrix = spot[7] * spotMatrix
                    
                    integrationMask, ringMask = imageSums_utilities.buildMasks(spotMatrix, 
                                                                               0, 
                                                                               ellipse_multiplicative_factor, 
                                                                               sigma_x, 
                                                                               sigma_y)
                            
                    ### CALCULATE var_ring_pix_img
                    n_ring = ringMask.sum()
                    integratedRing = numpy.multiply(ringMask, spotMatrix).sum()
                    avg_integratedRing = integratedRing / n_ring
                    
                    sum_sq = 0
                    n = 0                    
                    for i in range(0, spotMatrix.shape[0]):
                        for j in range(0, spotMatrix.shape[1]):
                            if ringMask[i, j] == 1:
                                n = n+1
                                bg_value = spotMatrix[i, j]
                                diff = bg_value - avg_integratedRing
                                sq = diff**2
                                sum_sq = sum_sq + sq  
                    if n != n_ring:
                        print 'PROBLEM'
                    var_ring_pix_img = sum_sq/n               
                    var_ring_pix_img_list.append(var_ring_pix_img)
                    
                                
        avg_var_ring_pix_img = numpy.average(var_ring_pix_img_list)
        
        flog.write('%5d %5d %5.2f %10d %7d %7d %9.2f %9.2f %d\n'%(h_label, k_label, q, n_terms, n_ellipse, n_ring, I_img_ellipse_ph, avg_var_ring_pix_img, n_bgPixels))
        
        # Use n_bgPixels (n pixels used in bg plane calculation) instead of n_ring
        r = n_ellipse/n_bgPixels
        noise =  numpy.sqrt( I_img_ellipse_ph + (r + 1) * n_ellipse * avg_var_ring_pix_img ) 
        signal_to_noise = numpy.sqrt(n_terms) * I_img_ellipse_ph / noise
        min_n_terms = (1/I_img_ellipse_ph) * ( 1 + (1/I_img_ellipse_ph) * n_ellipse * (r + 1) * avg_var_ring_pix_img )
        
        print signal_to_noise, min_n_terms
        print '\n'
                
        q_list.append(q)
        signal_to_noise_list.append(signal_to_noise)
        min_n_terms_list.append(min_n_terms)
        
#        # CONTROL PLOTS
#        controlPlots_flag = 0
#        if controlPlots_flag == 1:
#            myFigure = matplotlib.pyplot.figure()
#            myFigure.suptitle('%s  q = %.2f\n n_terms = %d n_ellipse = %d n_ring = %d\n I_img_ellipse_ph = %.2f\n avg_sigma_ring_pix_img = %.2f n_bgPixels = %d'%(label, q, 
#                                                                                                                                            n_terms, n_ellipse, n_ring, 
#                                                                                                                                            I_img_ellipse_ph, 
#                                                                                                                                            avg_var_ring_pix_img,
#                                                                                                                                            n_bgPixels), 
#                                                                                                                                            fontsize=6)
#            ax1 = myFigure.add_subplot(121)
#            ax1.tick_params(axis='both', which='major', labelsize=6)
#            ax1.set_title('MODULE DISPLACE - SUM', fontsize=5)
#            im1 = ax1.imshow(bgSubtractedOrbitSumMatrix_rotated, origin='lower', interpolation='nearest')     
#            divider1 = make_axes_locatable(ax1)    
#            cax1 = divider1.append_axes("right", size="10%", pad=0.05)  
#            cb1 = matplotlib.pyplot.colorbar(im1, cax=cax1)
#            for t in cb1.ax.get_yticklabels():
#                t.set_fontsize(6)    
#        
#            ax2 = myFigure.add_subplot(122)
#            ax2.tick_params(axis='both', which='major', labelsize=6)
#            ax2.set_title('MASKS', fontsize=5)
#            im2 = ax2.imshow(integrationMask+ringMask, origin='lower', interpolation='nearest') 
#            im2.set_cmap('Blues')   
#    
#            matplotlib.pyplot.tight_layout()
#            matplotlib.pyplot.savefig('%s/module_displace_%d_%d.png'%(outputFolder, h_label, k_label), dpi=4*96)
#            matplotlib.pyplot.close()

    qs_file = open('%s/signalToNoise%s_qs.pkl'%(outputFolder, folder_label), 'wb')
    pickle.dump(q_list, qs_file)
    qs_file.close()
    
    StoN_file = open('%s/signalToNoise%s_StoN.pkl'%(outputFolder, folder_label), 'wb')
    pickle.dump(signal_to_noise_list, StoN_file)
    StoN_file.close()
    
    minN_file = open('%s/signalToNoise%s_minN.pkl'%(outputFolder, folder_label), 'wb')
    pickle.dump(min_n_terms_list, minN_file)
    minN_file.close()

    flog.close()

if __name__ == "__main__":
    print "\n**** CALLING signalToNoise_module_displace ****"
    signalToNoise_module_displace(sys.argv[1:])