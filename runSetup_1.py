# -*- coding: utf-8 -*-
import os
import sys
import random

from runSetup_settings import tiltAngles
from runSetup_settings import nGoodFractions


# SETUP FOR CURRENT RUN
if len(sys.argv) != 2:
    print("[USAGE] %s runnumber" % sys.argv[0])
    sys.exit(-1)

runNumber = '%04d' % int(sys.argv[1])

# OUTPUT FOLDER
outFolder = './Output_r%s'%runNumber
if not os.path.exists(outFolder):
    os.mkdir(outFolder)



# CHECK CHEETAH PREPROCESSING RESULTS
# python checkH5content.py --runNumber <runNumber> --label <label>



clusterPath = '/das/work/p17/p17340/Cecilia_Casadei/2D-SFX'

# MAKE IMAGE LIST
imagesDirectoryName = '%s/UNIX_@_LCLS/r%s-images/data1'%(clusterPath, runNumber)

flag = 0
if flag == 1:
    os.system('python makeImageList.py --runNumber %s --imagesDirectoryName %s'
               %(runNumber, imagesDirectoryName))



# STORE IMAGE OBJECTS
tiltAngle = float(tiltAngles['%s'%runNumber])              # degrees 
imageListDirectory = '%s/ImageLists'%outFolder

flag = 0
if flag == 1:
    os.system('python storeImageObjects.py --runNumber %s \
                                           --tiltAngle %f \
                                           --imageListDirectory %s'
                                           %(runNumber, 
                                             tiltAngle, 
                                             imageListDirectory))



# Check image dictionary from command line: 
# python verify_pklContent.py --myPklFile <myPklFile> \
#                             --myImageFileName <myImageFileName>

# Display assembled images from command line:
# python verify_assembledImages.py --runNumber <runNumber> \
#                                  --imageNumber <imageNumber> \
#                                  --imageFolderName <imageFolderName>

# Browse assembled images from the command line:
# python verify_browseAssembledImages.py --runNumber <runNumber> \
#                                        --imageFolderName <imageFolderName>



# EXTRACT INFO FROM CHEETAH peaks.txt
selectedImageList = '%s/ImageLists/r%s_ImageNumbers_Filenames.txt'%(outFolder, 
                                                                    runNumber)
peaksFile = '%s/UNIX_\@_LCLS/r%s-good-modified-11/peaks.txt'%(clusterPath, 
                                                              runNumber)
geometryFile = '%s/Geometry/geometry.h5'%clusterPath # same for all runs
pixelSize = 0.000110  # m

flag = 0
if flag == 1:
    os.system('python extractExpInfo.py --runNumber %s \
                                        --selectedImageList %s \
                                        --peaksFile %s \
                                        --geometryFile %s \
                                        --pixelSize %f'
                                        %(runNumber, 
                                          selectedImageList, 
                                          peaksFile, 
                                          geometryFile, 
                                          pixelSize))



# Check image dictionary from command line: 
# python verify_pklContent.py --myPklFile <myPklFile> \
#                             --myImageFileName <myImageFileName>

# Display Cheetah peaks from command line:
# python verify_showCheetahPeaks.py --runNumber <runNumber> \
#                                   --imageNumber <imageNumber> \
#                                   --imageFolderName <imageFolderName> \
#                                   --geometryFile <geometryFile>



# BUILD REFERENCE RECIPROCAL LATTICE
referenceCellSize = 62.45    # A
hmax = 100                   # int
kmax = 100                   # int
resolutionLimit = 6.0        # A

flag = 0
if flag == 1:
    os.system('python reciprocalLattice.py --referenceCellSize %f \
                                           --hmax %d \
                                           --kmax %d \
                                           --resolutionLimit %f \
                                           --outFolder %s'
                                           %(referenceCellSize, 
                                             hmax, 
                                             kmax, 
                                             resolutionLimit, 
                                             outFolder))
    
    
    
# Check reciprocal lattice from command line:
# python verify_reciprocalLattice.py --reciprocalLatticeFile <reciprocalLatticeFile>
    
    
    
# INDEX LATTICES
detectorDistance      = 0.285  # m
radialTolerance       = 8      # pxls
pixelTolerance        = 12     # pxls
azimuthTolerance      = 3      # degrees
minNofPeaksPerLattice = 18     # int
maxNofPeaksPerImage   = 250    # int

indexingDirectory = '%s/LatticeIndexing'%outFolder

### INSTRUCTIONS TO USE MPI ###
# AT THE BEGINNING:
# salloc -N 5 -p day
# AT THE END:
# exit        
flag = 0
if flag == 1:   
    if not os.path.exists(indexingDirectory):   
            os.mkdir(indexingDirectory)      
    os.system('mpirun python latticeIndexing_mpi.py --referenceCellSize %f \
                                                    --runNumber %s \
                                                    --detectorDistance %f \
                                                    --pixelSize %f \
                                                    --radialTolerance %f \
                                                    --pixelTolerance %f \
                                                    --azimuthTolerance %f \
                                                    --minNofPeaksPerLattice %d \
                                                    --maxNofPeaksPerImage %d'
                                                    %(referenceCellSize, 
                                                      runNumber, 
                                                      detectorDistance, 
                                                      pixelSize, 
                                                      radialTolerance, 
                                                      pixelTolerance, 
                                                      azimuthTolerance, 
                                                      minNofPeaksPerLattice, 
                                                      maxNofPeaksPerImage))

    
 
flag = 0
if flag == 1:
    os.system('python latticeIndexing_mpi_merge.py --runNumber %s'%runNumber)
                                                          
                                                          
                                                            
flag = 0
if flag == 1:
    figuresPath = '%s/Figures'%indexingDirectory
    if not os.path.exists(figuresPath):
        os.mkdir(figuresPath)
    os.system('mpirun python latticeIndexing_mpi_plot.py --runNumber %s \
                                                          --detectorDistance %f \
                                                          --pixelSize %f'
                                                          %(runNumber, 
                                                            detectorDistance, 
                                                            pixelSize))
    os.system('rm %s/latticeDictionary*.pkl'%indexingDirectory)
       
    
# ORIENTATION AND CELL SIZE REFINEMENT
nSizeRefSteps            = 21
nOrientationRefSteps     = 21
widthSizeRefSteps        = 0.004
widthOrientationRefSteps = 0.2

refinementFolder = '%s/OrientationAndCellRefinement'%outFolder
    
flag = 0
if flag == 1:
       
    if not os.path.exists(refinementFolder):
        os.mkdir(refinementFolder)
    refinedPatternsFolder = '%s/RefinedPredictedPatterns'%refinementFolder
    if not os.path.exists(refinedPatternsFolder):
        os.mkdir(refinedPatternsFolder)
    figuresFolder = '%s/Figures'%refinementFolder
    if not os.path.exists(figuresFolder):
        os.mkdir(figuresFolder)
    
    os.system('mpirun python orientationAndCellRefinement_mpi.py --referenceCellSize %f \
                                                                 --runNumber %s \
                                                                 --nSizeRefSteps %d \
                                                                 --nOrientationRefSteps %d \
                                                                 --widthSizeRefSteps %f \
                                                                 --widthOrientationRefSteps %f \
                                                                 --resolutionLimit %f'
                                                                 %(referenceCellSize, 
                                                                   runNumber, 
                                                                   nSizeRefSteps, 
                                                                   nOrientationRefSteps, 
                                                                   widthSizeRefSteps, 
                                                                   widthOrientationRefSteps,
                                                                   resolutionLimit))
                                                                   
flag = 0
if flag == 1:
     os.system('python orientationAndCellRefinement_mpi_merge.py --runNumber %s'
                                                                  %runNumber)
               
               
               
# PLOT REFINED LATTICES
flag = 0
if flag == 1:            
    os.system('mpirun python plotRefinedLattices.py --runNumber %s'%runNumber)



# IMAGE PROCESSING
bgSubtractionMethod       = 'plane'
minimizationMethod        = '4Dbf'    # 4Dbf or Powell
lowResLimit               = 55.0
highResLimit              = 6.0
nCountsPerPhoton          = 26
integrationRadius         = 5
fractionDetectedThreshold = 0.45

flag = 0
if flag == 1:
    os.system('mpirun python processingAndIntegration.py --runNumber %s \
                                                         --bgSubtractionMethod %s \
                                                         --minimizationMethod %s \
                                                         --lowResLimit %f \
                                                         --highResLimit %f \
                                                         --nCountsPerPhoton %d \
                                                         --integrationRadius %d \
                                                         --geometryFile %s \
                                                         --imageFolder %s \
                                                         --fractionDetectedThreshold %f'
                                                         %(runNumber, 
                                                           bgSubtractionMethod, 
                                                           minimizationMethod, 
                                                           lowResLimit, 
                                                           highResLimit, 
                                                           nCountsPerPhoton, 
                                                           integrationRadius, 
                                                           geometryFile, 
                                                           imagesDirectoryName, 
                                                           fractionDetectedThreshold))
               
               
               
# Check individual lattice processing with single spot figures:
# python processingAndIntegration.py --runNumber <runNumber> --bgSubtractionMethod <bgSubtractionMethod>'
#                                    --minimizationMethod <minimizationMethod> --fractionDetectedThreshold <fractionDetectedThreshold> 
#                                    --lowResLimit <lowResLimit> --highResLimit <highResLimit>
#                                    --nCountsPerPhoton <nCountsPerPhoton> --integrationRadius <integrationRadius> --geometryFile <geometryFile>'
#                                    --imageFolder <imageFolder> --imageSelection <imageSelection> --latticeSelection <latticeSelection>' 
               
               
               
# MAKE LIST OF SPOT MATRICES
flag = 1
if flag == 1:
    os.system('python transform_makeSpotsMatrix.py --runNumber %s'%runNumber)
    
    
    
# DETERMINE TRANSFORMATIONS (spots to 7.1 A 2D-resolution are used)
deltaQrodThreshold = 0.005
n_minThreshold     = 6
nSeeds             = 20
nUsedLattices      = 'all'
nTriangles         = 100
nGoodFraction      = float(nGoodFractions['%s'%runNumber])   
    

flag = 0
if flag == 1:
    os.system('python transform_CCmethod_main.py --runNumber %s \
                                                 --dQrod %f \
                                                 --nMin %d \
                                                 --nSeeds %d \
                                                 --nLattices %s \
                                                 --nTriangles %d \
                                                 --nGoodFraction %f'
                                                 %(runNumber, 
                                                   deltaQrodThreshold, 
                                                   n_minThreshold, 
                                                   nSeeds, 
                                                   nUsedLattices, 
                                                   nTriangles, 
                                                   nGoodFraction))
                                                   

flag = 1
if flag == 1:
    for i in range(nSeeds):
        seed = random.sample(range(349), 1)
        os.system('mpirun python transform_CCmethod_mpi.py --runNumber %s \
                                                           --dQrod %f \
                                                           --nMin %d \
                                                           --nTriangles %d \
                                                           --nGoodFraction %f \
                                                           --seed %d'
                                                           %(runNumber, 
                                                             deltaQrodThreshold, 
                                                             n_minThreshold,
                                                             nTriangles, 
                                                             nGoodFraction,
                                                             seed[0]))
                                                             
        os.system('python transform_CCmethod_mpi_merge.py --runNumber %s \
                                                          --seed_id %d'
                                                          %(runNumber, 
                                                            i))
                                                            
    os.system('python transform_CCmethod_mpi_mergeSeeds.py --runNumber %s'
                                                           %(runNumber))
                        
                  
              
              
# DETERMINE TRANSFORMATIONS - SEEDS COMPARISON
flag = 1
if flag == 1:
    os.system('python transform_seedComparison.py --runNumber %s'%runNumber)
               

               
# APPLY TRANSFORMATIONS
flag = 0
if flag == 1:
    os.system('python transform_applyTransformations.py --runNumber %s'
                                                         %runNumber)
    
    
    
# SCALING
resolution_3D = 7.0 # A
n_minThreshold = 8

flag = 0
if flag == 1:
    os.system('python scaling.py --runNumber %s \
                                 --resolution_3D %f \
                                 --n_minThreshold %d'
                                 %(runNumber, 
                                   resolution_3D, 
                                   n_minThreshold))
    
    
    
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
    os.system('python plotRods.py --runNumber %s \
                                  --resolutionLimit %f'
                                  %(runNumber, 
                                    highResLimit))
