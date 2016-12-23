# -*- coding: utf-8 -*-
import joblib
import h5py
import numpy
import os
import pickle

import scipy.interpolate
import scipy.optimize

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot

import makeOrbits
import imageSums_utilities
import detectorModules

    

def imageSums():
    
    # PARAMETERS
    selectedRun = '0127'
    resolutionLimit = 4.0
    halfWidth = 15
    
    # FOLDERS
    outputFolder = './Output_imageSums_modules'    
    if not os.path.exists('%s'%outputFolder):
        os.mkdir('%s'%outputFolder)
    if not os.path.exists('%s/high_res_7_6'%outputFolder):
        os.mkdir('%s/high_res_7_6'%outputFolder)
    if not os.path.exists('%s/high_res_6_5'%outputFolder):
        os.mkdir('%s/high_res_6_5'%outputFolder)
    if not os.path.exists('%s/high_res_5_4'%outputFolder):
        os.mkdir('%s/high_res_5_4'%outputFolder)
    if not os.path.exists('%s/low_res'%outputFolder):
        os.mkdir('%s/low_res'%outputFolder)
            
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
    
    # EXTRACT DETECTOR MODULES BOUNDARIES
    modules = numpy.zeros((xGeometry_np.shape[0], xGeometry_np.shape[1]), dtype=numpy.int)
    modules = detectorModules.discontinuityMask(xGeometry_np, yGeometry_np, modules)
    
    row_interface = []
    row_interface.append(0)
    for i in range(0, modules.shape[0]):
        if modules[i, 0] == 1:
            row_interface.append(i)
    row_interface.append(modules.shape[0]-1)     #[0, 184, 369, 554, 739, 924, 1109, 1294, 1479]
    
    column_interface = []
    column_interface.append(0)
    for j in range(0, modules.shape[1]):
        if modules[0, j] == 1:
            column_interface.append(j)
    column_interface.append(modules.shape[1]-1)  #[0, 193, 387, 581, 775, 969, 1163, 1357, 1551]
    
    print row_interface
    print column_interface
    
    # RECIPROCAL CELL
    cellSize = 62.45
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I       
    
    # MAKE ORBIT OBJECTS LIST
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)                   # 220 orbits to 4 A
   
    # LOGGING
    fOpen = open('%s/log.txt'%outputFolder, 'a')
    fOpen.write('h, k, module bottom, module top, module left, module right, nTerms, Gauss_x0, Gauss_y0, sigma_x, sigma_y, refined_amplitude, Gauss_I\n')
    
    h_selected = 10
    k_selected = 3
   
    # FOR EACH ORBIT, FOR EACH MODULE, DO SECTOR SUMS
    for orbit in orbits:

        indices_1 = orbit.orbitIndices[0]
        indices_2 = orbit.orbitIndices[1]
        indices_3 = orbit.orbitIndices[2]
        h_1 = indices_1[0]
        k_1 = indices_1[1]
        h_2 = indices_2[0]
        k_2 = indices_2[1]
        h_3 = indices_3[0]
        k_3 = indices_3[1]
        
        # EXTRACT ORBIT LABEL:
        label = orbit.label
        h_label = label[0]
        k_label = label[1]
        
        if not (label == [h_selected, k_selected]):
            continue
        
        # CALCULATE q
        reciprocalVector = [h_label, k_label]*reciprocalCellRows
        q_x = reciprocalVector[0,0]         # A^(-1)
        q_y = reciprocalVector[0,1]         # A^(-1)
        q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
        
        if q < 0.8:
            continue
        
        print 'ORBIT: ', orbit.orbitIndices, 'LABEL: ', orbit.label
        
        results = [] # module top, module right, sigma_x, sigma_y, x0, y0, N, gauss_I
        nTot = 0
        
        # LOOP ON ALL MODULES
        for module_row in range(0, len(row_interface)-1):
            for module_column in range(0, len(column_interface)-1):
                bottom_bound = row_interface[module_row]
                top_bound = row_interface[module_row+1]
                left_bound = column_interface[module_column]
                right_bound = column_interface[module_column+1]
                print bottom_bound, top_bound, '---', left_bound, right_bound
                
                orbitSumMatrix_rotated = numpy.zeros(shape = (2*halfWidth, 2*halfWidth))              # 30x30
                orbitSumMatrix_rotated = numpy.matrix(orbitSumMatrix_rotated)
                
                orbitSumMatrix_nonRotated = numpy.zeros(shape = (2*halfWidth, 2*halfWidth))           # 30x30
                orbitSumMatrix_nonRotated = numpy.matrix(orbitSumMatrix_nonRotated)
                
                nTerms = 0
                   
                # LOOP ON ALL PROCESSED LATTICES 
                for index in range(0, len(lattices_names)):
                    lattice_name = lattices_names[index]
                    runNumber = lattice_name[10:14]
                    imageNumber = int(lattice_name[80:84])
                    latticeNumber = lattice_name[92]
                    latticeMatrix = lattices[index]                    # h_t k_t q_rod I_scaled flag i_unassembled j_unassembled scale
                    image_name = images_names[imageNumber-1].split()[1]
                    
                    if latticeMatrix[0, 4] == 0:                       # Check transformation and scaling flag
                        continue
                    
                    print runNumber, imageNumber, latticeNumber, image_name
                    
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
                            
                            if not (bottom_bound <= i < top_bound):
                                continue
                            if not (left_bound <= j < right_bound):
                                continue
                            
                            left_edge   = j - halfWidth
                            right_edge  = j + halfWidth
                            bottom_edge = i - halfWidth
                            top_edge    = i + halfWidth
                            if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                                continue
                            
                            ### EXTRACT 30x30 SECTOR ###
                            spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]
                            
                            ### SCALE ###
                            spotMatrix = numpy.asarray(spotMatrix, dtype=numpy.float32)
                            spotMatrix = spot[7] * spotMatrix
                            
                            ### DETERMINE AZIMUTH ON DETECTOR ###
                            xDetector = xGeometry_np[i, j]
                            yDetector = yGeometry_np[i, j]
                            if xDetector != 0:
                                detectorAzimuth = numpy.arctan(yDetector/xDetector)
                                if xDetector < 0 and yDetector > 0:
                                    detectorAzimuth = detectorAzimuth + numpy.pi
                                if xDetector < 0 and yDetector < 0:
                                    detectorAzimuth = detectorAzimuth - numpy.pi
                            else:
                                if yDetector > 0:
                                    detectorAzimuth = + numpy.pi /2
                                else:
                                    detectorAzimuth = - numpy.pi /2
                                    
                            ### DETERMINE MODULE ROTATION ANGLE ###
                            deltaXgeo = float(xGeometry_np[i, j+1] - xGeometry_np[i, j])
                            deltaYgeo = float(yGeometry_np[i, j+1] - yGeometry_np[i, j])
            
                            if deltaXgeo != 0:
                                moduleRotation = numpy.arctan(deltaYgeo/deltaXgeo) # From module to lab frame
                                if deltaXgeo < 0 and deltaYgeo > 0:
                                    moduleRotation = moduleRotation + numpy.pi
                                elif deltaXgeo < 0 and deltaYgeo < 0:
                                    moduleRotation = moduleRotation - numpy.pi
                            else:
                                if deltaYgeo > 0:
                                    moduleRotation = numpy.pi / 2
                                else:
                                    moduleRotation = -numpy.pi / 2
                             
                            ### ROTATION ANGLE ###
                            rotationAngle = - detectorAzimuth + moduleRotation
                            rotationAngle = - rotationAngle ### DUE TO CLOCKWISE ROTATION FUNCTION !!! 
                            
                            ### ROTATE ###  ### CLOCKWISE !!! ###
                            x_windows_pix = range(-spotMatrix.shape[1]/2, +spotMatrix.shape[1]/2)  # -18, -17, ..., +17
                            y_windows_pix = range(-spotMatrix.shape[0]/2, +spotMatrix.shape[0]/2)  # -18, -17, ..., +17
            
                            [X_windows_pix, Y_windows_pix] = numpy.meshgrid(x_windows_pix, y_windows_pix)
                            X_windows_pix = numpy.asarray(X_windows_pix, dtype=numpy.float32)
                            Y_windows_pix = numpy.asarray(Y_windows_pix, dtype=numpy.float32)
                            
                            X_windows_pix_rotated = numpy.cos(rotationAngle)*X_windows_pix - numpy.sin(rotationAngle)*Y_windows_pix
                            Y_windows_pix_rotated = numpy.sin(rotationAngle)*X_windows_pix + numpy.cos(rotationAngle)*Y_windows_pix
                            
                            f = scipy.interpolate.interp2d(x_windows_pix, y_windows_pix, spotMatrix, kind='linear')
                            spotMatrix_rotated = numpy.zeros(spotMatrix.shape)
                            
                            for columnIndex in range(0, spotMatrix_rotated.shape[1]):
                                for rowIndex in range(0, spotMatrix_rotated.shape[0]):
                                    rotated_x = X_windows_pix_rotated[rowIndex, columnIndex]
                                    rotated_y = Y_windows_pix_rotated[rowIndex, columnIndex]
                                    rotated_f = f(rotated_x, rotated_y)
                                    spotMatrix_rotated[rowIndex, columnIndex] = rotated_f                       
                            
                            ### SUM ###
                            orbitSumMatrix_rotated    = orbitSumMatrix_rotated    + spotMatrix_rotated
                            orbitSumMatrix_nonRotated = orbitSumMatrix_nonRotated + spotMatrix
                            
                            nTerms = nTerms + 1
                
                nTot = nTot + nTerms
                if nTerms == 0:
                    continue
                                                
                ### DIVIDE BY N OF TERMS IN THE SUM ###
                orbitSumMatrix_rotated    = orbitSumMatrix_rotated    / nTerms
                
                ### CALCULATE BACKGROUND ###
                background_rotated    = imageSums_utilities.calculateBackground(orbitSumMatrix_rotated,    outputFolder, 'rotated',     h_label, k_label)
                          
                ### BG SUBTRACTION ###
                bgSubtractedOrbitSumMatrix_rotated    = orbitSumMatrix_rotated    - background_rotated
                
                ### SINGLE PLOT ###
                myFigureObject = matplotlib.pyplot.figure()
                myAxesImageObject = matplotlib.pyplot.imshow(bgSubtractedOrbitSumMatrix_rotated, origin='lower', interpolation='nearest')
                matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f'%(h_label, k_label, orbit.resolution))
                myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                if orbit.resolution > 7.0:
                    matplotlib.pyplot.savefig('%s/low_res/bgsub_rotatetd_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96 )
                if 6.0 < orbit.resolution <= 7.0:
                    matplotlib.pyplot.savefig('%s/high_res_7_6/bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                if 5.0 < orbit.resolution <= 6.0:
                    matplotlib.pyplot.savefig('%s/high_res_6_5/bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                if 4.0 < orbit.resolution <= 5.0:
                    matplotlib.pyplot.savefig('%s/high_res_5_4/bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                matplotlib.pyplot.close()  
                
                ### TRY GAUSSIAN FIT OF BG SUBTRACTED, ROTATED SUM ###
                try:
                    n_x = bgSubtractedOrbitSumMatrix_rotated.shape[1]
                    n_y = bgSubtractedOrbitSumMatrix_rotated.shape[0]
                    x = numpy.linspace(0, n_x-1, n_x)
                    y = numpy.linspace(0, n_y-1, n_y)
                    x, y = numpy.meshgrid(x, y)                      # column_index, row_index
                    
                    data = bgSubtractedOrbitSumMatrix_rotated.ravel()           
                    data = data.T
                    data = numpy.asarray(data)
                    data = data.flatten()
                    
                    initial_x = float(n_x)/2
                    initial_y = float(n_y)/2
                    initial_guess = (numpy.amax(bgSubtractedOrbitSumMatrix_rotated), initial_x, initial_y, 2.0, 2.0)
                    
                    #param_bounds=([0,0,0,0,0],[numpy.inf,n_x,n_y, numpy.inf, numpy.inf])
                    popt, pcov = scipy.optimize.curve_fit(imageSums_utilities.twoD_Gaussian_simple, (x, y), data, p0=initial_guess) #bounds=param_bounds
                    data_fitted = imageSums_utilities.twoD_Gaussian_simple((x, y), *popt)       
                    
                    refined_amplitude = popt[0]
                    refined_x0 = popt[1]
                    refined_y0 = popt[2]
                    refined_sigma_x = popt[3]
                    refined_sigma_y = popt[4]
                          
                    ### ANALYTICAL GAUSSIAN INTEGRAL ###
                    gauss_integral = 2 * numpy.pi * refined_amplitude * refined_sigma_x * refined_sigma_y

                    ### PLOT GAUSS FIT ###
                    myFigureObject = matplotlib.pyplot.figure()
                    myAxesImageObject = matplotlib.pyplot.imshow(data.reshape(2*halfWidth, 2*halfWidth), origin='lower', interpolation='nearest')
                    matplotlib.pyplot.gca().contour(x, y, data_fitted.reshape(2*halfWidth, 2*halfWidth), 4, colors='w')
                    matplotlib.pyplot.title('Orbit: %d %d \nResolution: %.2f Gauss integral: %.1f counts'%(h_label, k_label, orbit.resolution, gauss_integral))
                    myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                    
                    if orbit.resolution > 7.0:
                        matplotlib.pyplot.savefig('%s/low_res/gauss_bgsub_rotatetd_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                    if 6.0 < orbit.resolution <= 7.0:
                        matplotlib.pyplot.savefig('%s/high_res_7_6/gauss_bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                    if 5.0 < orbit.resolution <= 6.0:
                        matplotlib.pyplot.savefig('%s/high_res_6_5/gauss_bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                    if 4.0 < orbit.resolution <= 5.0:
                        matplotlib.pyplot.savefig('%s/high_res_5_4/gauss_bgsub_rotated_orbit_%d_%d_n_%d_top_%d_right_%d.png'%(outputFolder, h_label, k_label, nTerms, top_bound, right_bound), dpi = 2*96)
                    matplotlib.pyplot.close()  
                    
                except:
                    print 'Gaussian fit not possible'
                    gauss_integral = numpy.nan
                    refined_x0 = numpy.nan
                    refined_y0 = numpy.nan
                    refined_sigma_x = numpy.nan
                    refined_sigma_y = numpy.nan
                 
                # LOGGING
                fOpen.write('%4d%4d %6d%6d%6d%6d %6d %10.3f%10.3f %8.2f%8.2f %10.1f %10.1f\n'
                %(h_label, k_label, bottom_bound, top_bound, left_bound, right_bound, nTerms, refined_x0, refined_y0, refined_sigma_x, refined_sigma_y, refined_amplitude, gauss_integral))
                if nTerms >= 50 and abs(gauss_integral) > 2*26 and refined_amplitude > 0:
                    module_results = [top_bound, right_bound, refined_sigma_x, refined_sigma_y, refined_x0, refined_y0, nTerms, gauss_integral]
                    results.append(module_results)
                    
    
    # SAVE RESULTS FOR SINGLE h, k
    results = numpy.asarray(results)   # module top, module right, sigma_x, sigma_y, x0, y0, N, gauss_I
    results_file = open('%s/module_results_%d_%d_totTerms_%d.pkl'%(outputFolder, h_selected, k_selected, nTot), 'wb')
    pickle.dump(results, results_file)
    results_file.close()    

    fOpen.close()

if __name__ == "__main__":
    print "\n**** CALLING imageSums ****"
    imageSums()