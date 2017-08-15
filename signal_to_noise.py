# -*- coding: utf-8 -*-
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy
import pickle
import joblib
import h5py

from mpl_toolkits.axes_grid1 import make_axes_locatable 

import imageSums_utilities

def buildMasks(sector, multiplicative_factor, sigma_x, sigma_y):
    # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
    integrationMask = numpy.zeros((sector.shape))   
    ringMask        = numpy.zeros((sector.shape)) 
    
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    x_axis = multiplicative_factor * sigma_x
    y_axis = multiplicative_factor * sigma_y
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
    
    integrationMask[numpy.where(distance < 1)] = 1
    ringMask[numpy.where(distance > 1.5)] = 1
    ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
    
    return integrationMask, ringMask


###################
def signalToNoise():
    
    ###        
    ### IMAGE SUM - MODULE DISPLACE
    ###

    
    # PARAMETERS
    selectedRun = '0127'
    nCountsPerPhoton = 26
    ellipse_multiplicative_factor = 2.5
    
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
    fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements/orbits.pkl'%selectedRun
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
    flog = open('%s/log.txt'%outputFolder, 'w')
    flog.write('#h, k, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_var_ring_pix_img, n_bgPixels\n')
    fLog = open('%s/log_StoN.txt'%outputFolder, 'w')
    fLog.write('#h, k, q, n_terms, StoN, n_min\n')
    
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
        # EXTRACT IMAGE SUM 
        # CONVERT TO PHOTON NUMBER 
        # KEEP NORMALIZATION ON N_TERMS
        ###
        bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtracted_total_sum # (30, 30) float64
        bgSubtractedOrbitSumMatrix_rotated = bgSubtractedOrbitSumMatrix_rotated / nCountsPerPhoton
        
        
        
        ###
        # BUILD MASKS
        ###
        sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
        integrationMask, ringMask = buildMasks(bgSubtractedOrbitSumMatrix_rotated, ellipse_multiplicative_factor, sigma_x, sigma_y)
        
        # n_ellipse, n_ring
        n_ellipse = integrationMask.sum()  
        n_ring    = ringMask.sum()         
        
        # n_bgPixels used in bg plane calculation
        n_bgPixels = orbit.n_bgPixels     
        
        # Total intensity on the ellipse, per term, expressed in photons
        I_img_ellipse_ph = numpy.multiply(integrationMask, bgSubtractedOrbitSumMatrix_rotated).sum()
        
        
        
        ###
        # CALCULATE THE AVERAGE OF STD DEVIATIONS OF RING INTENSITY FROM SINGLE IMAGES avg_sigma_ring_pix_img
        ###
        var_ring_pix_img_list = []
        n_max = 100
        halfWidth = int(bgSubtractedOrbitSumMatrix_rotated.shape[0]/2)
        
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
                    spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]  #(30, 30) float32
                    
                    # CONVERT SINGLE IMGS IN N OF PH
                    spotMatrix = spotMatrix / nCountsPerPhoton
                    
                    ### SCALE ###
                    spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                    spotMatrix = spot[7] * spotMatrix
                    
                    ### CALCULATE sigma_ring_pix_img
                    integratedRing = numpy.multiply(ringMask, spotMatrix).sum()
                    avg_integratedRing = integratedRing / n_ring
                    ring_intensities = numpy.multiply(ringMask, spotMatrix)                    
                    ring_intensities = ring_intensities.flatten().T # 900 elements                    
                
                    sum_sq = 0
                    n = 0
                    for ring_pxl in ring_intensities:
                        if ring_pxl != 0: # Out of ring mask (v. likely)                            
                            n = n+1
                            diff = ring_pxl - avg_integratedRing
                            sq = diff**2
                            sum_sq = sum_sq + sq
                        
                    var_ring_pix_img = sum_sq/n                   
                    var_ring_pix_img_list.append(var_ring_pix_img)
        
        avg_var_ring_pix_img = numpy.average(var_ring_pix_img_list)

        flog.write('%5d %5d %5.2f %10.1f %7.1f %7.1f %9.2f %9.2f %d\n'%(h_label, k_label, q, n_terms, n_ellipse, n_ring, I_img_ellipse_ph, avg_var_ring_pix_img, n_bgPixels))
        
        # Use n_bgPixels (n pixels used in bg plane calculation) instead of n_ring
        r = n_ellipse/n_bgPixels
        signal = numpy.sqrt(n_terms) * I_img_ellipse_ph
        noise = numpy.sqrt(I_img_ellipse_ph + (r + 1) * n_ellipse * avg_var_ring_pix_img )
        signal_to_noise = signal / noise
        min_n_terms = (1/I_img_ellipse_ph) + (1/(I_img_ellipse_ph**2)) * (r + 1) * n_ellipse * avg_var_ring_pix_img
        
        q_list.append(q)
        signal_to_noise_list.append(signal_to_noise)
        min_n_terms_list.append(min_n_terms)
        
        fLog.write('%5d %5d %5.2f %10.1f %7.1f\n'%(h_label, k_label, q, signal_to_noise, min_n_terms))

        #CONTROL PLOTS      
        controlPlots_flag = 0
        if controlPlots_flag == 1:               
            myFigure = matplotlib.pyplot.figure()
            myFigure.suptitle('%s  q = %.2f\n n_terms = %d n_ellipse = %d n_ring = %d\n I_img_ellipse_ph = %.2f\n avg_var_ring_pix_img = %.2f n_bgPixels = %d'%(label, q, 
                                                                                                                                            n_terms, n_ellipse, n_ring, 
                                                                                                                                            I_img_ellipse_ph, 
                                                                                                                                            avg_var_ring_pix_img,
                                                                                                                                            n_bgPixels), 
                                                                                                                                            fontsize=6)
            ax1 = myFigure.add_subplot(121)
            ax1.tick_params(axis='both', which='major', labelsize=6)
            ax1.set_title('MODULE DISPLACE - SUM', fontsize=5)
            im1 = ax1.imshow(bgSubtractedOrbitSumMatrix_rotated, origin='lower', interpolation='nearest')     
            divider1 = make_axes_locatable(ax1)    
            cax1 = divider1.append_axes("right", size="10%", pad=0.05)  
            cb1 = matplotlib.pyplot.colorbar(im1, cax=cax1)
            for t in cb1.ax.get_yticklabels():
                t.set_fontsize(6)    
        
            ax2 = myFigure.add_subplot(122)
            ax2.tick_params(axis='both', which='major', labelsize=6)
            ax2.set_title('MASKS', fontsize=5)
            im2 = ax2.imshow(integrationMask+ringMask, origin='lower', interpolation='nearest') 
            im2.set_cmap('Blues')   
    
            matplotlib.pyplot.tight_layout()
            matplotlib.pyplot.savefig('%s/module_displace_%d_%d.png'%(outputFolder, h_label, k_label), dpi=4*96)
            matplotlib.pyplot.close()

    qs_file = open('%s/signalToNoise_qs.pkl'%outputFolder, 'wb')
    pickle.dump(q_list, qs_file)
    qs_file.close()
    
    StoN_file = open('%s/signalToNoise_StoN.pkl'%outputFolder, 'wb')
    pickle.dump(signal_to_noise_list, StoN_file)
    StoN_file.close()
    
    minN_file = open('%s/signalToNoise_minN.pkl'%outputFolder, 'wb')
    pickle.dump(min_n_terms_list, minN_file)
    minN_file.close()
    
    flog.close()
    fLog.close()

if __name__ == "__main__":
    print "\n**** CALLING signalToNoise ****"
    signalToNoise()