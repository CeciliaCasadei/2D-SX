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

def buildMasks(sector, precision_factor, multiplicative_factor, sigma_x, sigma_y):
    # PREPARE INTEGRATION MASK (ELLIPSE) AND RING MASK (ELLIPTICAL RING)
    integrationMask = numpy.zeros((sector.shape))   
    ringMask        = numpy.zeros((sector.shape)) 
    
    colIdx, rowIdx = numpy.meshgrid(numpy.arange(integrationMask.shape[1]), numpy.arange(integrationMask.shape[0]))
    x_axis = 10**precision_factor * multiplicative_factor * sigma_x
    y_axis = 10**precision_factor * multiplicative_factor * sigma_y
    centerX = integrationMask.shape[1] / 2
    centerY = integrationMask.shape[0] / 2
    distance = ((rowIdx-centerY)**2)/(y_axis**2) + ((colIdx-centerX)**2)/(x_axis**2)
    
    integrationMask[numpy.where(distance < 1)] = 1
    ringMask[numpy.where(distance > 1.5)] = 1
    ringMask[numpy.where(distance > 3.5)] = 0 #was 2.1
    
    return integrationMask, ringMask

def signalToNoise_no_module_displace():
    
    ###        
    ### IMAGE SUM - NO MODULE DISPLACE
    ###

    
    # PARAMETERS
    selectedRun = '0127'
    nCountsPerPhoton = 26
    ellipse_multiplicative_factor = 2.5
    precision_factor = 0
    k_squared = (10**precision_factor)**2
    
    # FOLDERS
    outputFolder = './Output_r%s/Output_signal_to_noise_no_module_displace'%selectedRun
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)       
    
    # EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
    fRead = open('./Output_r%s/Output_imageSums_sigmaFits/sigmaXCurveParameters.pkl'%selectedRun, 'rb')
    sigmaXCurveParameters = pickle.load(fRead)                                               
    fRead.close()
    
    fRead = open('./Output_r%s/Output_imageSums_sigmaFits/sigmaYCurveParameters.pkl'%selectedRun, 'rb')
    sigmaYCurveParameters = pickle.load(fRead)                                               
    fRead.close()   
    
    # LOAD ORBITS
    fileToOpen = './Output_r%s/Output_imageSums/orbits.pkl'%selectedRun
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
    
    # EXTRACT GEOMETRY
    geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5'     # same for all runs
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    
    # LOG
    flog = open('%s/log.txt'%outputFolder, 'w')
    flog.write('#h, k, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_sigma_ring_pix_img\n')
        
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
        # REMOVE NORMALIZATION ON N_TERMS
        ###
        bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtractedOrbitSumMatrix_rotated # (30, 30) float64
        bgSubtractedOrbitSumMatrix_rotated = bgSubtractedOrbitSumMatrix_rotated / nCountsPerPhoton
        bgSubtractedOrbitSumMatrix_rotated = n_terms * bgSubtractedOrbitSumMatrix_rotated
        
        
        ###
        # BUILD MASKS
        ###
        sigma_x = imageSums_utilities.quadratic(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.line_plus_sigmoid(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2], sigmaYCurveParameters[3], sigmaYCurveParameters[4])    
        integrationMask, ringMask = buildMasks(bgSubtractedOrbitSumMatrix_rotated, precision_factor, ellipse_multiplicative_factor, sigma_x, sigma_y)
        
        # n_ellipse, n_ring
        n_ellipse = integrationMask.sum()
        n_ring    = ringMask.sum()
        
        # Total intensity on the ellipse, in the sum, expressed in photons
        I_sum_ellipse_ph = numpy.multiply(integrationMask, bgSubtractedOrbitSumMatrix_rotated).sum()
        
        # Total intensity on the ellipse, per image, expressed in photons
        I_img_ellipse_ph = I_sum_ellipse_ph / n_terms
        
        ###
        # CALCULATE THE AVERAGE OF STD DEVIATIONS OF RING INTENSITY FROM SINGLE IMAGES avg_sigma_ring_pix_img
        ###
        halfWidth = bgSubtractedOrbitSumMatrix_rotated.shape[0] / 2   #15 pxls
        sigma_ring_pix_img_list = []
        n_max = 100
        
        # Loop on lattices
        for index in range(0, len(lattices_names)):
            
            if len(sigma_ring_pix_img_list) > n_max:
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
                    
                    ### DETERMINE AZIMUTH ON DETECTOR ###
                    detectorAzimuth = imageSums_utilities.calculate_detectorAzimuth(xGeometry_np, yGeometry_np, i, j)
                            
                    ### DETERMINE MODULE ROTATION ANGLE ###
                    moduleRotation = imageSums_utilities.calculate_moduleRotation(xGeometry_np, yGeometry_np, i, j)
                     
                    ### ROTATION ANGLE ###
                    rotationAngle = - detectorAzimuth + moduleRotation
                    rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                    
                    ### ROTATE ###  ### CLOCKWISE !!! ###
                    spotMatrix_rotated = imageSums_utilities.clockWiseRotation(spotMatrix, rotationAngle)                    
                    
                    ### CALCULATE sigma_ring_pix_img
                    integratedRing = numpy.multiply(ringMask, spotMatrix_rotated).sum()
                    avg_integratedRing = integratedRing / n_ring
                    ring_intensities = numpy.multiply(ringMask, spotMatrix_rotated)                    
                    ring_intensities = ring_intensities.flatten().T # 900 elements
                    
                
                    sum_sq = 0
                    n = 0
                    for ring_pxl in ring_intensities:
                        if ring_pxl != 0: # Out of ring mask (v. likely)                            
                            n = n+1
                            diff = ring_pxl - avg_integratedRing
                            sq = diff**2
                            sum_sq = sum_sq + sq
                    if n_ring != n:
                        print 'n_ring (%d) different from n (%d)'%(n_ring, n)
                    
                        
                    sigma_ring_pix_img = numpy.sqrt(sum_sq/n)                    
                    sigma_ring_pix_img_list.append(sigma_ring_pix_img)
        
        avg_sigma_ring_pix_img = numpy.average(sigma_ring_pix_img_list)
       
        print label, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_sigma_ring_pix_img
        print len(sigma_ring_pix_img_list)  
        flog.write('%5d %5d %5.2f %10.1f %7.1f %7.1f %10.2f %9.2f %9.2f\n'%(h_label, k_label, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_sigma_ring_pix_img))
        
        
        signal_to_noise = I_sum_ellipse_ph / numpy.sqrt( k_squared*( I_sum_ellipse_ph + n_terms*((1/n_ring)+(1/n_ellipse))*(n_ellipse**2)*(avg_sigma_ring_pix_img**2)  ) )
        min_n_terms = (k_squared/(I_img_ellipse_ph**2))*(I_img_ellipse_ph + ((1/n_ring)+(1/n_ellipse)) * (n_ellipse**2) * (avg_sigma_ring_pix_img**2) )
        
        q_list.append(q)
        signal_to_noise_list.append(signal_to_noise)
        min_n_terms_list.append(min_n_terms)

        # CONTROL PLOTS
#        myFigure = matplotlib.pyplot.figure()
#        myFigure.suptitle('%s  q = %.2f\n n_terms = %d n_ellipse = %d n_ring = %d\n I_sum_ellipse_ph = %.2f I_img_ellipse_ph = %.2f avg_sigma_ring_pix_img = %.2f'%(label, q, 
#                                                                                                                                        n_terms, n_ellipse, n_ring, 
#                                                                                                                                        I_sum_ellipse_ph, 
#                                                                                                                                        I_img_ellipse_ph, 
#                                                                                                                                        avg_sigma_ring_pix_img), fontsize=7)
#        ax1 = myFigure.add_subplot(121)
#        ax1.tick_params(axis='both', which='major', labelsize=6)
#        ax1.set_title('NO DISPLACE - SUM', fontsize=5)
#        im1 = ax1.imshow(bgSubtractedOrbitSumMatrix_rotated, origin='lower', interpolation='nearest')     
#        divider1 = make_axes_locatable(ax1)    
#        cax1 = divider1.append_axes("right", size="10%", pad=0.05)  
#        cb1 = matplotlib.pyplot.colorbar(im1, cax=cax1)
#        for t in cb1.ax.get_yticklabels():
#            t.set_fontsize(6)    
#    
#        ax2 = myFigure.add_subplot(122)
#        ax2.tick_params(axis='both', which='major', labelsize=6)
#        ax2.set_title('MASKS', fontsize=5)
#        im2 = ax2.imshow(integrationMask+ringMask, origin='lower', interpolation='nearest') 
#        im2.set_cmap('Blues')   
#
#        matplotlib.pyplot.tight_layout()
#        matplotlib.pyplot.savefig('%s/no_displace_%d_%d.png'%(outputFolder, h_label, k_label), dpi=4*96)
#        matplotlib.pyplot.close()


# BINNING    
    q_plot = []
    signal_to_noise_plot = [] 
    min_n_terms_plot = []

    bins = numpy.linspace(min(q_list), max(q_list), 25)
    for i in range(0, len(bins)-1):
        left_q  = bins[i]
        right_q = bins[i+1]
        q_bin               = [q_list[i]               for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]    
        signal_to_noise_bin = [signal_to_noise_list[i] for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
        min_n_terms_bin     = [min_n_terms_list[i]     for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
        
        print len(q_bin)
        q_avg_bin = numpy.average(q_bin)
        signal_to_noise_avg_bin = numpy.average(signal_to_noise_bin)
        min_n_terms_avg_bin = numpy.average(min_n_terms_bin)

        q_plot.append(q_avg_bin)
        signal_to_noise_plot.append(signal_to_noise_avg_bin)
        min_n_terms_plot.append(min_n_terms_avg_bin)

    matplotlib.pyplot.figure()     
    matplotlib.pyplot.scatter(q_list, signal_to_noise_list, s=3)
    matplotlib.pyplot.scatter(q_plot, signal_to_noise_plot, s=130, facecolors='r', edgecolors='r', alpha = 0.25)
    matplotlib.pyplot.gca().set_ylim([0.01, 1000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"Signal/Noise", fontsize = 12, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/signal_to_noise_no_displace.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    matplotlib.pyplot.figure()     
    matplotlib.pyplot.scatter(q_list, min_n_terms_list, s=3)
    matplotlib.pyplot.scatter(q_plot, min_n_terms_plot, s=130, facecolors='r', edgecolors='r', alpha = 0.25)
    matplotlib.pyplot.gca().set_ylim([0.001, 1000000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"Minimum number of terms", fontsize = 12, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/min_n_terms_no_displace.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    
    flog.close()
   
###############################
###        
### IMAGE SUM - MODULE DISPLACE
###
def signalToNoise_module_displace():
    import makeOrbits
    
    # PARAMETERS
    selectedRun = '0127'
    resolution_limit = 4
    nCountsPerPhoton = 26
    ellipse_multiplicative_factor = 2.5
    precision_factor = 1
    k_squared = (10**precision_factor)**2 #100
    halfWidth = 25
    truncated_halfWidth = 15
    
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
    
    # PREPARE LATTICES TO IMAGES MATCHING
    imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%selectedRun    
    lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                           %(selectedRun, selectedRun, selectedRun)) # 629 lattices  
    lattices_names = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(selectedRun, selectedRun, selectedRun), 'r')
    lattices_names = list(lattices_names)    
    images_names = open('./Output_r%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(selectedRun, selectedRun))
    images_names = list(images_names)
    
    # EXTRACT GEOMETRY
    geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5'     # same for all runs
    geometryData = h5py.File(geometryFile, 'r')
    xGeometry = geometryData['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryData['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
        
    ###
    # FOR EACH ORBIT, CALCULATE THE AVERAGE OF STD DEVIATIONS OF RING INTENSITY FROM SINGLE IMAGES avg_sigma_ring_pix_img
    ###    
    orbit_avg_sigma_ring_pix_img = []
    orbits = makeOrbits.makeOrbitsFunction(resolution_limit)
    for orbit in orbits:
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
        
        ###
        # BUILD MASKS
        ###
        sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
        sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
        integrationMask, ringMask = buildMasks(numpy.zeros(shape=(300, 300)), precision_factor, ellipse_multiplicative_factor, sigma_x, sigma_y)
        n_ring = ringMask.sum()
        
        sigma_ring_pix_img_list = []
        n_max = 100
        
        # Loop on lattices
        for index in range(0, len(lattices_names)):
            
            if len(sigma_ring_pix_img_list) > n_max:
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
                    
                    ### DETERMINE AZIMUTH ON DETECTOR ###
                    detectorAzimuth = imageSums_utilities.calculate_detectorAzimuth(xGeometry_np, yGeometry_np, i, j)
                            
                    ### DETERMINE MODULE ROTATION ANGLE ###
                    moduleRotation = imageSums_utilities.calculate_moduleRotation(xGeometry_np, yGeometry_np, i, j)
                     
                    ### ROTATION ANGLE ###
                    rotationAngle = - detectorAzimuth + moduleRotation
                    rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                    
                    ### ROTATE ###  ### CLOCKWISE !!! ###
                    spotMatrix_rotated = imageSums_utilities.clockWiseRotation(spotMatrix, rotationAngle)          
                    
                    x0 = float(halfWidth)
                    y0 = float(halfWidth)
                    expanded_spotMatrix_rotated = imageSums_utilities.recenter(spotMatrix_rotated, x0, y0, precision_factor, truncated_halfWidth) # 50x50 -> 300x300, NO BG SUB, NO NORMALIZATION
                    
                    ### CALCULATE sigma_ring_pix_img
                    integratedRing = numpy.multiply(ringMask, expanded_spotMatrix_rotated).sum()
                    avg_integratedRing = integratedRing / n_ring
                    ring_intensities = numpy.multiply(ringMask, expanded_spotMatrix_rotated)                    
                    ring_intensities = ring_intensities.flatten().T # 900 elements                    
                
                    sum_sq = 0
                    n = 0
                    for ring_pxl in ring_intensities:
                        if ring_pxl != 0: # Out of ring mask (v. likely)                            
                            n = n+1
                            diff = ring_pxl - avg_integratedRing
                            sq = diff**2
                            sum_sq = sum_sq + sq
                    if n_ring != n:
                        print 'n_ring (%d) different from n (%d)'%(n_ring, n)                    
                        
                    sigma_ring_pix_img = numpy.sqrt(sum_sq/n)                    
                    sigma_ring_pix_img_list.append(sigma_ring_pix_img)
        
        avg_sigma_ring_pix_img = numpy.average(sigma_ring_pix_img_list)
        print h_label, k_label, avg_sigma_ring_pix_img
        print len(sigma_ring_pix_img_list)  
        orbit_avg_sigma_ring_pix_img.append([h_label, k_label, avg_sigma_ring_pix_img])
    orbit_avg_sigma_ring_pix_img = numpy.asarray(orbit_avg_sigma_ring_pix_img)    
    
    #######################################################################################
    colors = ['b', 'm', 'c']
    orbit_folders = ['', '_100_lattices', '_10_lattices']
    matplotlib.pyplot.figure()
    
    StoN = []
    index = 0
    for suffix in orbit_folders:
        
        # LOAD ORBITS
        fileToOpen = './Output_r%s/Output_imageSums_moduleDisplacements%s/orbits.pkl'%(selectedRun, suffix)
        fRead = open(fileToOpen, 'rb')
        orbitsDictionary = pickle.load(fRead)                                               
        fRead.close()
        
        # LOG
        flog = open('%s/log_%s.txt'%(outputFolder, suffix), 'w')
        flog.write('#h, k, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_sigma_ring_pix_img, n_bgPixels\n')
            
        q_list = []
        signal_to_noise_list = []
        min_n_terms_list = []
        
        # LOOP ON ORBITS
        for key, orbit in orbitsDictionary.items():        
            label = orbit.label        
            
            # h, k
            h_label = label[0]
            k_label = label[1]
                        
            # q
            q = 2 * numpy.pi / orbit.resolution
            
            # n_terms
            n_terms = orbit.nTerms
           
            ###
            # EXTRACT IMAGE SUM 
            # CONVERT TO PHOTON NUMBER 
            # REMOVE NORMALIZATION ON N_TERMS
            ###
            bgSubtractedOrbitSumMatrix_rotated = orbit.bgSubtracted_total_sum # (300, 300) float64
            bgSubtractedOrbitSumMatrix_rotated = bgSubtractedOrbitSumMatrix_rotated / nCountsPerPhoton
            bgSubtractedOrbitSumMatrix_rotated = n_terms * bgSubtractedOrbitSumMatrix_rotated            
            
            ###
            # BUILD MASKS
            ###
            sigma_x = imageSums_utilities.poly_4(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
            sigma_y = imageSums_utilities.poly_4(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2])
            integrationMask, ringMask = buildMasks(bgSubtractedOrbitSumMatrix_rotated, precision_factor, ellipse_multiplicative_factor, sigma_x, sigma_y)
            
            # n_ellipse, n_ring
            n_ellipse = integrationMask.sum()
            n_ring    = ringMask.sum()
            
            # n_bgPixels used in bg plane calculation
            n_bgPixels = orbit.n_bgPixels
            
            # Total intensity on the ellipse, in the sum, expressed in photons
            I_sum_ellipse_ph = numpy.multiply(integrationMask, bgSubtractedOrbitSumMatrix_rotated).sum()
            
            # Total intensity on the ellipse, per image, expressed in photons
            I_img_ellipse_ph = I_sum_ellipse_ph / n_terms

            ###
            # CALCULATE THE AVERAGE OF STD DEVIATIONS OF RING INTENSITY FROM SINGLE IMAGES avg_sigma_ring_pix_img
            ###
            avg_sigma_ring_pix_img = numpy.nan
            for i in orbit_avg_sigma_ring_pix_img:
                if i[0] == h_label and i[1] == k_label:
                    avg_sigma_ring_pix_img = i[2]
                    break
            if numpy.isnan(avg_sigma_ring_pix_img):
                print 'PROBLEM'
                    
            
           
            print label, q, n_terms, n_ellipse, n_ring, I_sum_ellipse_ph, I_img_ellipse_ph, avg_sigma_ring_pix_img, n_bgPixels
            flog.write('%5d %5d %5.2f %10.1f %7.1f %7.1f %10.2f %9.2f %9.2f %d\n'%(h_label, k_label, q, n_terms, n_ellipse, n_ring, 
                                                                                   I_sum_ellipse_ph, I_img_ellipse_ph, 
                                                                                   avg_sigma_ring_pix_img, n_bgPixels))
            
            # Use n_bgPixels (n pixels used in bg plane calculation) instead of n_ring
            signal_to_noise = I_sum_ellipse_ph / numpy.sqrt( k_squared*( I_sum_ellipse_ph + n_terms*((1/n_bgPixels)+(1/n_ellipse))*(n_ellipse**2)*(avg_sigma_ring_pix_img**2)  ) )
            min_n_terms = (k_squared/(I_img_ellipse_ph**2))*(I_img_ellipse_ph + ((1/n_bgPixels)+(1/n_ellipse)) * (n_ellipse**2) * (avg_sigma_ring_pix_img**2) )
            
            q_list.append(q)
            signal_to_noise_list.append(signal_to_noise)
            min_n_terms_list.append(min_n_terms)
            
        # BINNING    
        q_plot = []
        signal_to_noise_plot = [] 
        min_n_terms_plot = []
    
        bins = numpy.linspace(min(q_list), max(q_list), 25)
        for i in range(0, len(bins)-1):
            left_q  = bins[i]
            right_q = bins[i+1]
            q_bin               = [q_list[i]               for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]    
            signal_to_noise_bin = [signal_to_noise_list[i] for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
            min_n_terms_bin     = [min_n_terms_list[i]     for i in range(0, len(q_list)) if left_q <= q_list[i] <= right_q]
            
            q_avg_bin = numpy.average(q_bin)
            signal_to_noise_avg_bin = numpy.average(signal_to_noise_bin)
            min_n_terms_avg_bin = numpy.average(min_n_terms_bin)
    
            q_plot.append(q_avg_bin)
            signal_to_noise_plot.append(signal_to_noise_avg_bin)
            min_n_terms_plot.append(min_n_terms_avg_bin)
            
       
        matplotlib.pyplot.scatter(q_plot, signal_to_noise_plot, marker='o', s=6, color=colors[index])
        #matplotlib.pyplot.scatter(q_plot, min_n_terms_plot, marker='o', s=6, color=colors[index])
        flog.close()
        index = index + 1
        StoN.append(signal_to_noise_plot)
        
    matplotlib.pyplot.gca().set_ylim([0.001, 1000])
    matplotlib.pyplot.gca().set_yscale('log')    
    matplotlib.pyplot.gca().set_xlabel(r"q ($\AA^{-1}$)", fontsize = 12, rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r"Signal/Noise", fontsize = 12, rotation = 'vertical')
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_compare.png'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.savefig('%s/signal_to_noise_module_displace_compare.pdf'%(outputFolder), dpi=4*96)
    #matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_compare.png'%(outputFolder), dpi=4*96)
    #matplotlib.pyplot.savefig('%s/min_n_terms_module_displace_compare.pdf'%(outputFolder), dpi=4*96)
    matplotlib.pyplot.close()
    
    StoN_all = StoN[0]
    StoN_100 = StoN[1]
    StoN_10  = StoN[2]
    
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
    
    ratios_all_to_100 = numpy.asarray(ratios_all_to_100).flatten()
    ratios_all_to_10 = numpy.asarray(ratios_all_to_10).flatten()
    ratios_100_to_10 = numpy.asarray(ratios_100_to_10).flatten()
    print 'Ratio 586 to 10', numpy.average(ratios_all_to_10)
    print 'Ratio 100 to 10', numpy.average(ratios_100_to_10)
    print 'Ratio 586 to 100', numpy.average(ratios_all_to_100)

if __name__ == "__main__":
    print "\n**** CALLING signalToNoise_module_displace ****"
    signalToNoise_module_displace()