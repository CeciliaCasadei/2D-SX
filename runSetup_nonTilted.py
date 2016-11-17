# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_nonTilted.py

# SETUP FOR CURRENT RUN

runNumber = '0127'

# OUTPUT FOLDER
outFolder = './Output_r%s'%runNumber
if not os.path.exists(outFolder):
    os.mkdir(outFolder)



# LOG CURRENT SETTINGS
os.system('cp ./runSetup_r0127.py ./Output_r%s/runSetup.log'%runNumber)



# CHECK CHEETAH PREPROCESSING RESULTS
# python checkH5content.py --runNumber <runNumber> --label <label>


# MAKE IMAGE LIST
imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%runNumber

flag = 0
if flag == 1:
    os.system('python makeImageList.py --runNumber %s --imagesDirectoryName %s'
               %(runNumber, imagesDirectoryName))



# STORE IMAGE OBJECTS
tiltAngle          = 0              # degrees 
imageListDirectory = '%s/ImageLists'%outFolder

flag = 0
if flag == 1:
    os.system('python storeImageObjects.py --runNumber %s --tiltAngle %f --imageListDirectory %s'
               %(runNumber, tiltAngle, imageListDirectory))



# Check image dictionary from command line: 
# python verify_pklContent.py --myPklFile <myPklFile> --myImageFileName <myImageFileName>

# Display assembled images from command line:
# python verify_assembledImages.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName>

# Browse assembled images from the command line:
# python verify_browseAssembledImages.py --runNumber <runNumber> --imageFolderName <imageFolderName>



# EXTRACT INFO FROM CHEETAH peaks.txt
selectedImageList = '%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(outFolder, runNumber)
peaksFile         = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_\@_LCLS/r%s-good-modified-9/peaks.txt'%runNumber
geometryFile      = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5' # same for all runs
pixelSize         = 0.000110 # m

flag = 0
if flag == 1:
    os.system('python extractExpInfo.py --runNumber %s --selectedImageList %s --peaksFile %s --geometryFile %s --pixelSize %f'
               %(runNumber, selectedImageList, peaksFile, geometryFile, pixelSize))



# Check image dictionary from command line: 
# python verify_pklContent.py --myPklFile <myPklFile> --myImageFileName <myImageFileName>

# Display Cheetah peaks from command line:
# python verify_showCheetahPeaks.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName> --geometryFile <geometryFile>



# BUILD REFERENCE RECIPROCAL LATTICE
referenceCellSize = 62.45    # A
hmax = 100                   # int
kmax = 100                   # int
resolutionLimit = 4.0        # A

flag = 0
if flag == 1:
    os.system('python reciprocalLattice.py --referenceCellSize %f --hmax %d --kmax %d --resolutionLimit %f --outFolder %s'
               %(referenceCellSize, hmax, kmax, resolutionLimit, outFolder))
    
    
    
# Check reciprocal lattice from command line:
# python verify_reciprocalLattice.py --reciprocalLatticeFile <reciprocalLatticeFile>
    
    
    
# INDEX LATTICES
detectorDistance      = 0.235     # m
radialTolerance       = 6         # pxls
pixelTolerance        = 10        # pxls
azimuthTolerance      = 3         # degrees
minNofPeaksPerLattice = 20        # int
maxNofPeaksPerImage   = 250       # int
        
flag = 0
if flag == 1:
    os.system('python latticeIndexing.py --referenceCellSize %f --runNumber %s --detectorDistance %f --pixelSize %f --radialTolerance %f --pixelTolerance %f --azimuthTolerance %f --minNofPeaksPerLattice %d --maxNofPeaksPerImage %d --geometryFile %s'
    %(referenceCellSize, runNumber, detectorDistance, pixelSize, radialTolerance, pixelTolerance, azimuthTolerance, minNofPeaksPerLattice, maxNofPeaksPerImage, geometryFile))
    
    
    
# ORIENTATION AND CELL SIZE REFINEMENT
nSizeRefSteps            = 25
nOrientationRefSteps     = 25
widthSizeRefSteps        = 0.002
widthOrientationRefSteps = 0.2
    
flag = 0
if flag == 1:
    os.system('python orientationAndCellRefinement.py --referenceCellSize %f --runNumber %s --nSizeRefSteps %d --nOrientationRefSteps %d --widthSizeRefSteps %f --widthOrientationRefSteps %f --hmax %d --kmax %d --resolutionLimit %f'
               %(referenceCellSize, runNumber, nSizeRefSteps, nOrientationRefSteps, widthSizeRefSteps, widthOrientationRefSteps, hmax, kmax, resolutionLimit))
               
               
               
# PLOT REFINED LATTICES
flag = 0
if flag == 1:
    os.system('python plotRefinedLattices.py --runNumber %s'%runNumber)
    
    
    
# PLOT REFINED LATTICES - PAPER FIGURE
flag = 0
if flag == 1:
    os.system('python plotRefinedLattices_imageOverlap.py --runNumber %s'%runNumber)



# IMAGE PROCESSING
bgSubtractionMethod       = 'plane'
minimizationMethod        = '4Dbf'    # 4Dbf or Powell
lowResLimit               = 55.0      # A
highResLimit              = 4.0       # A
nCountsPerPhoton          = 26
integrationRadius         = 5         # pxls
fractionDetectedThreshold = 0.28

flag = 0
if flag == 1:
    os.system('python processingAndIntegration.py --runNumber %s --bgSubtractionMethod %s --minimizationMethod %s --lowResLimit %f --highResLimit %f --nCountsPerPhoton %d --integrationRadius %d --geometryFile %s --imageFolder %s --fractionDetectedThreshold %f'
               %(runNumber, bgSubtractionMethod, minimizationMethod, lowResLimit, highResLimit, nCountsPerPhoton, integrationRadius, geometryFile, imagesDirectoryName, fractionDetectedThreshold))
               
               
               
# Check individual lattice processing with single spot figures:
# python processingAndIntegration.py --runNumber <runNumber> --bgSubtractionMethod <bgSubtractionMethod>'
#                                    --minimizationMethod <minimizationMethod> --fractionDetectedThreshold <fractionDetectedThreshold> 
#                                    --lowResLimit <lowResLimit> --highResLimit <highResLimit>
#                                    --nCountsPerPhoton <nCountsPerPhoton> --integrationRadius <integrationRadius> --geometryFile <geometryFile>'
#                                    --imageFolder <imageFolder> --imageSelection <imageSelection> --latticeSelection <latticeSelection>' 


               
# PLOT CELL SIZE DISTRIBUTION 
flag = 0
if flag == 1:
    os.system('python cellSizeDistribution.py --runNumber %s'%runNumber)
    

               
# MAKE LIST OF SPOT MATRICES
flag = 0
if flag == 1:
    os.system('python transform_makeSpotsMatrix.py --runNumber %s'%runNumber)
    
    
    
# DETERMINE TRANSFORMATIONS
deltaQrodThreshold = 0.001
n_minThreshold = 6
nSeeds = 6
nUsedLattices = 'all'
nTriangles = 100
nGoodFraction = 0.7
    
flag = 0
if flag == 1:
    os.system('python transform_CCmethod_main.py --runNumber %s --dQrod %f --nMin %d --nSeeds %d --nLattices %s --nTriangles %d --nGoodFraction %f'
              %(runNumber, deltaQrodThreshold, n_minThreshold, nSeeds, nUsedLattices, nTriangles, nGoodFraction))
              
              
              
# DETERMINE TRANSFORMATIONS - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python transform_seedComparison.py --runNumber %s --nSeeds %d --dQrod %f'
               %(runNumber, nSeeds, deltaQrodThreshold))
               

               
# APPLY TRANSFORMATIONS
flag = 0
if flag == 1:
    os.system('python transform_applyTransformations.py --runNumber %s'%runNumber)
    
    
    
# SCALING
flag = 0
if flag == 1:
    os.system('python scaling.py --runNumber %s --dQrod %f'%(runNumber, deltaQrodThreshold))
    
    
    
# SCALING - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python scaling_seedComparison.py --runNumber %s'%runNumber)
    


# PLOT SCALE FACTOR DISTRIBUTION 
flag = 0
if flag == 1:
    os.system('python scaleFactorDistribution.py --runNumber %s'%runNumber)    


    
# SCALING - APPLY SCALES
flag = 0
if flag == 1:
    os.system('python scaling_applyScales.py --runNumber %s'%runNumber)
        
    
    
# PLOT RODS AND ORBIT HISTOGRAMS
flag = 0
if flag == 1:
    os.system('python plotRods_r0127.py --runNumber %s'%runNumber)
        
    
    
# PLOT PAPER FIGURE (4 MERGING HISTOGRAMS)
flag = 0
if flag == 1:
    os.system('python plotSelectedMergingHistograms.py --runNumber %s'%runNumber)



# PLOT MERGING EFFICIENCY
flag = 0
if flag == 1:
    os.system('python mergingEfficiency.py --runNumber %s'%runNumber)
    


### ************************************ IMAGE SUMS ************************************ ###
    
# IMAGE SUMS
flag = 0
if flag == 1:
    os.system('python imageSums.py')

    
    
####### ALTERNATIVE #######
# IMAGE SUMS - FIRST BG SUBTRACTION AND THEN SUM
#flag = 0
#if flag == 1:
#    os.system('python imageSums_bgSubAndSum.py')
    
    

# IMAGE SUMS - PAPER PLOTS
flag = 0
if flag == 1:
    os.system('python imageSums_plots.py')
    
    
    
# IMAGE SUMS - SIGMAS vs Q PLOTS
flag = 0
if flag == 1:
    os.system('python imageSums_sigmaVsQ.py')
    
    
    
# IMAGE SUMS - PAPER PLOTS
flag = 0
if flag == 1:
    os.system('python imageSums_sigmaVsQ_paperFigure.py')
    
    
    
# IMAGE SUMS - PAPER PLOTS
flag = 0
if flag == 1:
    os.system('python imageSums_ellipse_Gauss_comparison.py')
    
    
    
# ELLIPTICAL INTEGRATION ON SINGLE IMAGES
flag = 0
if flag == 1:
    os.system('python singleImage_ellipse_integration.py')
    

    
# SCALE ELLIPTICAL INTEGRATED IMAGES
ellipticalIntegrationFolder = './Output_elliptical_integration'
flag = 0
if flag == 1:
    os.system('python scaling.py --runNumber %s --dQrod %f --outputFolder %s'%(runNumber, deltaQrodThreshold, ellipticalIntegrationFolder))
        


# SCALING - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python scaling_seedComparison.py --runNumber %s --outputFolder %s'%(runNumber, ellipticalIntegrationFolder))
        
        
        
# SCALING - APPLY SCALES
flag = 0
if flag == 1:
    os.system('python scaling_applyScales.py --runNumber %s --outputFolder %s'%(runNumber, ellipticalIntegrationFolder))
    
    
    
# PLOT RODS AND ORBIT HISTOGRAMS
flag = 0
if flag == 1:
    os.system('python plotRods_r0127.py --runNumber %s --outputFolder %s'%(runNumber, ellipticalIntegrationFolder))
    
    
    
# PLOT PAPER FIGURE (4 MERGING HISTOGRAMS)
flag = 0
if flag == 1:
    os.system('python plotSelectedMergingHistograms.py --runNumber %s --outputFolder %s'%(runNumber, ellipticalIntegrationFolder))
    


# COMPARE METHODS: MERGING AFTER ELLIPTICAL INTEGRATION OF SINGLE IMAGES - GAUSS INTEGRAL OF IMAGE SUMS    
flag = 0
if flag == 1:
    os.system('python compare_I_no_sum_ellipse_I_sum_gauss.py')
    
    
    
# CALCULATE RESOLUTION CUTOFF USING SINGLE IMAGE METHOD
flag = 0
if flag == 1:
    os.system('python resolutionCutoff_singleImage_ellipseIntegration.py')
    
    
    
# PLOT RESOLUTION CUTOFF USING SINGLE IMAGE METHOD AND IMAGE SUMS METHOD
flag = 1
if flag == 1:
    os.system('python resolutionCutoff_plots.py')