# -*- coding: utf-8 -*-
import os

# INSTRUCTIONS: 
# python runSetup_nonTilted.py

# SETUP FOR CURRENT RUN

runNumber = '0052'

# OUTPUT FOLDER
outFolder = './Output_r%s'%runNumber
if not os.path.exists(outFolder):
    os.mkdir(outFolder)



# LOG CURRENT SETTINGS
os.system('cp ./runSetup_nonTilted_anomalous.py ./Output_r%s/runSetup_nonTilted_anomalous.log'%runNumber)



# CHECK CHEETAH PREPROCESSING RESULTS
# python checkH5content.py --runNumber <runNumber> --label <label>


# MAKE IMAGE LIST
imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%runNumber
#imagesDirectoryName = '/mnt/das-gpfs/home/casadei_c/work/casadei/UNIX_@_LCLS/r%s-images/data1'%runNumber

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
peaksFile         = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_\@_LCLS/r%s-good-modified-11/peaks.txt'%runNumber
#peaksFile = '/mnt/das-gpfs/home/casadei_c/work/casadei/UNIX_\@_LCLS/r%s-good-modified-11/peaks.txt'%runNumber
geometryFile      = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5' # same for all runs
#geometryFile = '/mnt/das-gpfs/home/casadei_c/work/casadei/Geometry/geometry.h5' # same for all runs
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
detectorDistance      = 0.351     # m
radialTolerance       = 6         # pxls
pixelTolerance        = 10        # pxls
azimuthTolerance      = 3         # degrees
minNofPeaksPerLattice = 20        # int
maxNofPeaksPerImage   = 300       # int
        
flag = 0
if flag == 1:
    os.system('python latticeIndexing.py --referenceCellSize %f --runNumber %s --detectorDistance %f --pixelSize %f --radialTolerance %f --pixelTolerance %f --azimuthTolerance %f --minNofPeaksPerLattice %d --maxNofPeaksPerImage %d'
    %(referenceCellSize, runNumber, detectorDistance, pixelSize, radialTolerance, pixelTolerance, azimuthTolerance, minNofPeaksPerLattice, maxNofPeaksPerImage))
    
    
    
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
highResLimit              = 5.8       # A
nCountsPerPhoton          = 26
integrationRadius         = 5         # pxls
fractionDetectedThreshold = 0.40

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
productThreshold = 0.45
flag = 0
if flag == 1:
    os.system('python scaling.py --runNumber %s --dQrod %f --productThreshold %f'%(runNumber, deltaQrodThreshold, productThreshold))
    
    
    
# SCALING - SEEDS COMPARISON
flag = 0
if flag == 1:
    os.system('python scaling_seedComparison.py --runNumber %s'%runNumber)
    


# PLOT SCALE FACTOR DISTRIBUTION 
flag = 0
if flag == 1:
    os.system('python scaleFactorDistribution.py --runNumber %s'%runNumber)    


    
# SCALING - APPLY SCALES
flag = 1
if flag == 1:
    os.system('python scaling_applyScales.py --runNumber %s'%runNumber)
        
    
    
# PLOT RODS AND ORBIT HISTOGRAMS
flag = 0
if flag == 1:
    os.system('python plotRods_zero_tilt.py --runNumber %s'%runNumber)
        
    
    
# PLOT PAPER FIGURE (4 MERGING HISTOGRAMS)
flag = 0
if flag == 1:
    os.system('python plotSelectedMergingHistograms.py --runNumber %s'%runNumber)

flag = 0
if flag == 1:
    os.system('python plotSelectedMergingHistograms_scalingEffect.py --runNumber %s'%runNumber)

# PLOT MERGING EFFICIENCY
flag = 0
if flag == 1:
    os.system('python mergingEfficiency.py --runNumber %s'%runNumber)
    


### ************************************ IMAGE SUMS ************************************ ###
### ****************************** NO MODULE TRANSLATION ******************************* ###
    
# IMAGE SUMS
flag = 0
if flag == 1:
    os.system('python imageSums.py')

    
    
####### ALTERNATIVE #######
#IMAGE SUMS - FIRST BG SUBTRACTION AND THEN SUM
#flag = 0
#if flag == 1:
#    os.system('python imageSums_bgSubAndSum.py')
### REALIZE THAT MANY GAUSS FITS ARE BAD -> THIS WAY OF BG SUBTRACTING IS NOT ADEQUATE
    

# IMAGE SUMS - PAPER PLOTS
flag = 0
if flag == 1:
    os.system('python imageSums_singleOrbitPlots.py')
    
    
    
# IMAGE SUMS - SIGMAS vs Q BEHAVIOUR
flag = 0
if flag == 1:
    os.system('python imageSums_sigmaVsQ.py')
### REALIZE THAT THERE IS A PROBLEM IN CALIBRATION OF MODULE POSITIONS -> NEED TO INTRO TRANSLATION CORRECTIONS
    
    
    
### ************************************ IMAGE SUMS ************************************ ###
### ******************************** MODULE TRANSLATION ******************************** ###
flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements.py')
    
    
flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements_extract.py')
    

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules.py')
        

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_sigmaVsQ.py')
    
    

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_finalIntegration.py')
    # prepare_ellipseIntegralLatexTable.py    

    
    
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_integrationPlots.py')



### Data quality evaluation ###   
N_lattices = 10
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_N_lattices.py --N_lattices %d'%N_lattices)


    
N_lattices = 100
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_N_lattices.py --N_lattices %d'%N_lattices)
    
 
   
flag = 0
if flag == 1:
    os.system('python signal_to_noise.py')
    
    
    
flag = 0
if flag == 1:
    os.system('python signal_to_noise_plot.py')
    
    

flag = 0
if flag == 1:
    os.system('python signal_to_noise_compare.py')



N_lattices = 10
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py --N_lattices %d'%N_lattices)
  
  
    
N_lattices = 100
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py --N_lattices %d'%N_lattices)
    
    

N_lattices = 586
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py --N_lattices %d'%N_lattices)


 
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_calculateCC.py')
    
    
    
# PAPER FIGURE WITH IMG SUMS
flag = 0
if flag == 1:
    os.system('python imageSums_displacedModules_plots.py')