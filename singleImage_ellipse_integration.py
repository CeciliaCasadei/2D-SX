# -*- coding: utf-8 -*-
import pickle
import numpy
import os
import joblib
import h5py
import scipy.interpolate

import imageSums_utilities

# EXTRACT SIGMA_X AND SIGMA_Y CURVE PARAMETERS
fRead = open('./Output_imageSums_sigmaFits/sigmaXCurveParameters.pkl', 'rb')
sigmaXCurveParameters = pickle.load(fRead)                                               
fRead.close()

fRead = open('./Output_imageSums_sigmaFits/sigmaYCurveParameters.pkl', 'rb')
sigmaYCurveParameters = pickle.load(fRead)                                               
fRead.close()

# PARAMETERS
selectedRun = '0127'
halfWidth = 15
ellipse_multiplicative_factor = 4.0 #was 2.5

# FOLDERS
outputFolder = './Output_elliptical_integration'    
if not os.path.exists('%s'%outputFolder):
    os.mkdir('%s'%outputFolder)
        
# PREPARE LATTICES TO IMAGES MATCHING
imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%selectedRun    
lattices = joblib.load('./Output_r%s/transformAndScale/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'
                       %(selectedRun, selectedRun, selectedRun)) # 586 lattices  
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

print 'LATTICE NAMES LIST LENGTH: %d'%len(lattices_names)
print 'LATTICES LIST LENGTH: %d'%len(lattices)

# RECIPROCAL CELL
cellSize = 62.45
directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
reciprocalCellRows = 2 * numpy.pi * directCell.I  

# LOOP ON ALL PROCESSED LATTICES 
ellipseIntegratedLattices_list = []

for index in range(0, len(lattices_names)):
    ellipseIntegratedLattice = []
    
    lattice_name = lattices_names[index]
    runNumber = lattice_name[10:14]
    imageNumber = int(lattice_name[80:84])
    latticeNumber = lattice_name[92]
    latticeMatrix = lattices[index]                    # h_t k_t q_rod I_unscaled flag i_unassembled j_unassembled 
    image_name = images_names[imageNumber-1].split()[1]
    print runNumber, imageNumber, latticeNumber, image_name
       
    # LOAD UNASSEMBLED IMAGE
    unassembledDataFile = h5py.File('%s/%s'%(imagesDirectoryName, image_name), 'r')
    unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
    unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)         #### !!!!! ####  (1480, 1552)
    nRows = unassembledData.shape[0]
    nColumns = unassembledData.shape[1]
    
    # IN THE CURRENT LATTICE, LOOP ON SPOTS 
    for spot in latticeMatrix:
        ellipseIntegratedSpot = []
        h = int(spot[0])
        k = int(spot[1])
        if spot[4] == 0: # flag
            ellipseIntegratedSpot = [spot[0], spot[1], spot[2], spot[3], spot[4], spot[5], spot[6]]
        else:
            reciprocalVector = [h, k]*reciprocalCellRows
            q_x = reciprocalVector[0,0]         # A^(-1)
            q_y = reciprocalVector[0,1]         # A^(-1)
            q = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
            
            sigma_x = imageSums_utilities.quadratic(q, sigmaXCurveParameters[0], sigmaXCurveParameters[1], sigmaXCurveParameters[2])
            sigma_y = imageSums_utilities.line_plus_sigmoid(q, sigmaYCurveParameters[0], sigmaYCurveParameters[1], sigmaYCurveParameters[2], sigmaYCurveParameters[3], sigmaYCurveParameters[4])
            
            i = int(spot[5])
            j = int(spot[6])
            left_edge   = j - halfWidth
            right_edge  = j + halfWidth
            bottom_edge = i - halfWidth
            top_edge    = i + halfWidth
            if left_edge < 0 or right_edge > nColumns or bottom_edge < 0 or top_edge > nRows:
                ellipseIntegratedSpot = [spot[0], spot[1], spot[2], numpy.nan, spot[4], spot[5], spot[6]]
            else:
                ### EXTRACT 30x30 SECTOR ###
                spotMatrix = unassembledData[bottom_edge:top_edge, left_edge:right_edge]
                
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
                
                # BACKGROUND SUBTRACTION
                background_rotated  = imageSums_utilities.calculateBackground_noImg(spotMatrix_rotated)
                bgSubtracted_spotMatrix_rotated = spotMatrix_rotated - background_rotated
                
                # SUM ON ELLIPSE MASK
                ellipse_intensity = imageSums_utilities.integrate_ellipse(bgSubtracted_spotMatrix_rotated, sigma_x, sigma_y, ellipse_multiplicative_factor)
                
                ellipseIntegratedSpot = [spot[0], spot[1], spot[2], ellipse_intensity, spot[4], spot[5], spot[6]]
        ellipseIntegratedLattice.append(ellipseIntegratedSpot)
    ellipseIntegratedLattice = numpy.asarray(ellipseIntegratedLattice, dtype=numpy.float32)
    ellipseIntegratedLattices_list.append(ellipseIntegratedLattice)
    
print 'LIST LENGTH: %d'%len(ellipseIntegratedLattices_list)

if not os.path.exists('%s/spotsMatricesList-Transformed-r%s'%(outputFolder, selectedRun)):
    os.mkdir('%s/spotsMatricesList-Transformed-r%s'%(outputFolder, selectedRun))
joblib.dump(ellipseIntegratedLattices_list, '%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, selectedRun, selectedRun))