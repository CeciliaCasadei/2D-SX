# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_1_rXXXX.py
# python runSetup_1_rXXXX.py > log.log



# SETUP FOR CURRENT RUN

runNumber = '0200'

# OUTPUT FOLDER
outFolder = './Output_r%s'%runNumber
if not os.path.exists(outFolder):
    os.mkdir(outFolder)



# LOG CURRENT SETTINGS
os.system('cp ./runSetup_1_r0200.py ./Output_r%s/runSetup.log'%runNumber)



# CHECK CHEETAH PREPROCESSING RESULTS
# python checkH5content.py --runNumber <runNumber> --label <label>


# MAKE IMAGE LIST
imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%runNumber

flag = 1
if flag == 1:
    os.system('python makeImageList.py --runNumber %s --imagesDirectoryName %s'
               %(runNumber, imagesDirectoryName))



# STORE IMAGE OBJECTS
tiltAngle = 20               # degrees 
imageListDirectory = '%s/ImageLists'%outFolder

flag = 1
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
peaksFile = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_\@_LCLS/r%s-good-modified-11/peaks.txt'%runNumber
geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5' # same for all runs
pixelSize = 0.000110         # m

flag = 1
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
resolutionLimit = 5.0        # A

flag = 1
if flag == 1:
    os.system('python reciprocalLattice.py --referenceCellSize %f --hmax %d --kmax %d --resolutionLimit %f --outFolder %s'
               %(referenceCellSize, hmax, kmax, resolutionLimit, outFolder))
    
    
    
# Check reciprocal lattice from command line:
# python verify_reciprocalLattice.py --reciprocalLatticeFile <reciprocalLatticeFile>
    
    
    
# INDEX LATTICES
detectorDistance = 0.285     # m
radialTolerance = 8          # pxls
pixelTolerance = 12          # pxls
azimuthTolerance = 3         # degrees
minNofPeaksPerLattice = 18   # int
maxNofPeaksPerImage = 250    # int
        
flag = 1
if flag == 1:
    os.system('python latticeIndexing.py --referenceCellSize %f --runNumber %s --detectorDistance %f --pixelSize %f --radialTolerance %f --pixelTolerance %f --azimuthTolerance %f --minNofPeaksPerLattice %d --maxNofPeaksPerImage %d --geometryFile %s'
    %(referenceCellSize, runNumber, detectorDistance, pixelSize, radialTolerance, pixelTolerance, azimuthTolerance, minNofPeaksPerLattice, maxNofPeaksPerImage, geometryFile))
    
    
    
# ORIENTATION AND CELL SIZE REFINEMENT
nSizeRefSteps = 21
nOrientationRefSteps = 21
widthSizeRefSteps = 0.004
widthOrientationRefSteps = 0.2
    
flag = 1
if flag == 1:
    os.system('python orientationAndCellRefinement.py --referenceCellSize %f --runNumber %s --nSizeRefSteps %d --nOrientationRefSteps %d --widthSizeRefSteps %f --widthOrientationRefSteps %f --hmax %d --kmax %d --resolutionLimit %f'
               %(referenceCellSize, runNumber, nSizeRefSteps, nOrientationRefSteps, widthSizeRefSteps, widthOrientationRefSteps, hmax, kmax, resolutionLimit))
               
               
               
# PLOT REFINED LATTICES
flag = 0
if flag == 1:
    os.system('python plotRefinedLattices.py --runNumber %s'%runNumber)



# IMAGE PROCESSING
bgSubtractionMethod = 'plane'
minimizationMethod = '4Dbf'    # 4Dbf or Powell
lowResLimit = 55.0
highResLimit = 7.1
nCountsPerPhoton = 26
integrationRadius = 5
fractionDetectedThreshold = 0.55

flag = 1
if flag == 1:
    os.system('python processingAndIntegration.py --runNumber %s --bgSubtractionMethod %s --minimizationMethod %s --lowResLimit %f --highResLimit %f --nCountsPerPhoton %d --integrationRadius %d --geometryFile %s --imageFolder %s --fractionDetectedThreshold %f'
               %(runNumber, bgSubtractionMethod, minimizationMethod, lowResLimit, highResLimit, nCountsPerPhoton, integrationRadius, geometryFile, imagesDirectoryName, fractionDetectedThreshold))
               
               
               
# Check individual lattice processing with single spot figures:
# python processingAndIntegration.py --runNumber <runNumber> --bgSubtractionMethod <bgSubtractionMethod>'
#                                    --minimizationMethod <minimizationMethod> --fractionDetectedThreshold <fractionDetectedThreshold> 
#                                    --lowResLimit <lowResLimit> --highResLimit <highResLimit>
#                                    --nCountsPerPhoton <nCountsPerPhoton> --integrationRadius <integrationRadius> --geometryFile <geometryFile>'
#                                    --imageFolder <imageFolder> --imageSelection <imageSelection> --latticeSelection <latticeSelection>' 
               
               
               
# MAKE LIST OF SPOT MATRICES
flag = 0
if flag == 1:
    os.system('python transform_makeSpotsMatrix.py --runNumber %s'%runNumber)
    
    
    
# DETERMINE TRANSFORMATIONS
deltaQrodThreshold = 0.005
deltaSelfThreshold = 0
n_minThreshold = 6
nSeeds = 6
nUsedLattices = 'all'
nTriangles = 100
nGoodFraction = 0.6
    
flag = 0
if flag == 1:
    os.system('python transform_CCmethod_main.py --runNumber %s --dQrod %f --dSelf %f --nMin %d --nSeeds %d --nLattices %s --nTriangles %d --nGoodFraction %f'
              %(runNumber, deltaQrodThreshold, deltaSelfThreshold, n_minThreshold, nSeeds, nUsedLattices, nTriangles, nGoodFraction))
              
              
              
# DETERMINE TRANSFORMATIONS - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python transform_seedComparison.py --runNumber %s --nSeeds %d --dQrod %f --dSelf %f'
               %(runNumber, nSeeds, deltaQrodThreshold, deltaSelfThreshold))
               

               
# APPLY TRANSFORMATIONS
flag = 0
if flag == 1:
    os.system('python transform_applyTransformations.py --runNumber %s'%runNumber)
    
    
    
# SCALING
flag = 0
if flag == 1:
    os.system('python scaling.py --runNumber %s'%runNumber)
    
    
    
# SCALING - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python scaling_seedComparison.py --runNumber %s'%runNumber)
    
    
    
# SCALING - APPLY SCALES
flag = 0
if flag == 1:
    os.system('python scaling_applyScales.py --runNumber %s'%runNumber)
        
    
    
# PLOT RODS
flag = 0
if flag == 1:
    os.system('python plotRods.py --runNumber %s'%runNumber)