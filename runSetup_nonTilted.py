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
os.system('cp ./runSetup_nonTilted.py ./Output_r%s/runSetup_nonTilted.log'%runNumber)



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
peaksFile         = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_\@_LCLS/r%s-good-modified-9/peaks.txt'%runNumber
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
detectorDistance      = 0.235     # m
radialTolerance       = 6         # pxls
pixelTolerance        = 10        # pxls
azimuthTolerance      = 3         # degrees
minNofPeaksPerLattice = 20        # int
maxNofPeaksPerImage   = 250       # int
        
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
### *********************** MODULE TRANSLATION BY INTERPOLATION ************************ ###
    
halfwidth = 15
approximateWavelength = 1.48 #A

flag = 0
if flag == 1:
    os.system('python calculate_partialSums.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d         \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f'
              %(runNumber, highResLimit, 
                halfwidth, imagesDirectoryName, 
                geometryFile, referenceCellSize, hmax, kmax, 
                tiltAngle, approximateWavelength, detectorDistance, pixelSize))



moduleDisplace_I_threshold = 2          # photons
moduleDisplace_nTerms_threshold = 50
moduleDisplace_sigma_threshold = 3.5    # pxls
moduleDisplace_distance_threshold = 7.5 # pxls

flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements.py \
               --selectedRun %s --resolutionLimit %f \
               --intensity_threshold %f --nTerms_threshold %d --sigma_threshold %f --distance_threshold %f \
               --nCountsPerPhoton %d --geometryFile %s'
               %(runNumber, highResLimit,
                 moduleDisplace_I_threshold, moduleDisplace_nTerms_threshold, moduleDisplace_sigma_threshold, moduleDisplace_distance_threshold,
                 nCountsPerPhoton, geometryFile))

    
    
flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements_extract.py --selectedRun %s --halfWidth %d --geometryFile %s'
              %(runNumber, halfwidth, geometryFile))


    
lowFluctuationThreshold = 1.5 # STANDARD IS 2.0

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules.py --selectedRun %s --resolutionLimit %f --halfWidth %d --lowFluctuationThreshold %f'
              %(runNumber, highResLimit, halfwidth, lowFluctuationThreshold))
        


sigmaVsQ_phThreshold = 1.0
sigmaVsQ_dThreshold  = 10.0

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_sigmaVsQ.py --selectedRun %s --photonThreshold %f --dThreshold %f --nCountsPerPhoton %d --cellSize %f'
              %(runNumber, sigmaVsQ_phThreshold, sigmaVsQ_dThreshold, nCountsPerPhoton, referenceCellSize))
    
    

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_finalIntegration.py')
    


flag = 0
if flag == 1:
    os.system('python prepare_ellipseIntegralLatexTable.py --selectedRun %s'%runNumber)



flag = 0  
if flag == 1:
    os.system('python imageSums_displaced_modules_plots.py')
    
    

flag = 0
if flag == 1:
    os.system('python plot_I_ellipse_vs_q.py --selectedRun %s --cellSize %f'%(runNumber, referenceCellSize))


    
flag = 0
if flag == 1:
    os.system('python signal_to_noise.py')
    
    
    
flag = 0
if flag == 1:
    os.system('python signal_to_noise_plot.py')
    
    


    
### ************************************ IMAGE SUMS ************************************ ###
### ******************** MODULE TRANSLATION BY PIXEL CONVERSION ************************ ###

halfwidth = 25
approximateWavelength = 1.48 #A

flag = 0
if flag == 1:
    os.system('python calculate_partialSums.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d         \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f'
              %(runNumber, highResLimit, 
                halfwidth, imagesDirectoryName, 
                geometryFile, referenceCellSize, hmax, kmax, 
                tiltAngle, approximateWavelength, detectorDistance, pixelSize))



moduleDisplace_I_threshold = 2          # photons
moduleDisplace_nTerms_threshold = 50
moduleDisplace_sigma_threshold = 4      # pxls
moduleDisplace_distance_threshold = 15  # pxls

flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements_total_I_threshold.py \
               --selectedRun %s --resolutionLimit %f \
               --intensity_threshold %f --nTerms_threshold %d --sigma_threshold %f --distance_threshold %f \
               --nCountsPerPhoton %d --geometryFile %s'
               %(runNumber, highResLimit,
                 moduleDisplace_I_threshold, moduleDisplace_nTerms_threshold, moduleDisplace_sigma_threshold, moduleDisplace_distance_threshold,
                 nCountsPerPhoton, geometryFile))



flag = 0
if flag == 1:
    os.system('python calculate_moduleDisplacements_extract.py --selectedRun %s --halfWidth %d --geometryFile %s'
              %(runNumber, halfwidth, geometryFile))



lowFluctuationThreshold = 2.0 # STANDARD IS 2.0
precisionFactor = 1

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_pixelConversion.py --selectedRun %s --resolutionLimit %f --halfWidth %d --lowFluctuationThreshold %f --precisionFactor %d'
              %(runNumber, highResLimit, halfwidth, lowFluctuationThreshold, precisionFactor))


              
sigmaVsQ_phThreshold = 1.0
sigmaVsQ_dThreshold  = 1000000000

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_sigmaVsQ.py --selectedRun %s --photonThreshold %f --dThreshold %f --nCountsPerPhoton %d --cellSize %f'
              %(runNumber, sigmaVsQ_phThreshold, sigmaVsQ_dThreshold, nCountsPerPhoton, referenceCellSize))



ellipse_multiplicative_factor = 2.5

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_finalIntegration_pixelConversion.py --selectedRun %s --nCountsPerPhoton %d --ellipse_multiplicative_factor %f --precisionFactor %d'
              %(runNumber, nCountsPerPhoton, ellipse_multiplicative_factor, precisionFactor))      



flag = 0
if flag == 1:
    os.system('python prepare_ellipseIntegralLatexTable.py --selectedRun %s'%runNumber)
    
    
    
flag = 0
if flag == 1:
    os.system('python plot_I_ellipse_vs_q.py --selectedRun %s --cellSize %f'%(runNumber, referenceCellSize))
    
    
    
flag = 0 
if flag == 1:
    os.system('python imageSums_displaced_modules_plots_pixelConversion.py --selectedRun %s --ellipse_multiplicative_factor %f --precisionFactor %d'
              %(runNumber, ellipse_multiplicative_factor, precisionFactor))



d_threshold = 3.0 # pxls

flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_integrationPlots.py --selectedRun %s --d_threshold %f --halfWidth %d'
              %(runNumber, d_threshold, halfwidth))    
    


label = "!"    
flag = 0
if flag == 1:
    os.system('python signal_to_noise_pixelConversion.py \
              --selectedRun %s --nCountsPerPhoton %d \
              --ellipse_multiplicative_factor %f --precisionFactor %d \
              --halfWidth %d --label %s'
              %(runNumber, nCountsPerPhoton, 
                ellipse_multiplicative_factor, precisionFactor, 
                halfwidth, label))




nBins = 25

flag = 0
if flag == 1:
    os.system('python signal_to_noise_plot_pixelConversion.py --selectedRun %s --nBins %d'
    %(runNumber, nBins))


    
### SIGNAL TO NOISE vs Q, WITH DIFFERENT N LATTICES ###   
    
N_lattices = 10
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_N_lattices.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d         \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f \
               --precision_factor %d --N_lattices %d'
               %(runNumber, highResLimit,
                 halfwidth, imagesDirectoryName,
                 geometryFile, referenceCellSize, hmax, kmax,
                 tiltAngle, approximateWavelength, detectorDistance, pixelSize,
                 precisionFactor, N_lattices))



N_lattices = 100
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_N_lattices.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d         \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f \
               --precision_factor %d --N_lattices %d'
               %(runNumber, highResLimit,
                 halfwidth, imagesDirectoryName,
                 geometryFile, referenceCellSize, hmax, kmax,
                 tiltAngle, approximateWavelength, detectorDistance, pixelSize,
                 precisionFactor, N_lattices))
    
 

label = "_10_lattices"    
flag = 0
if flag == 1:
    os.system('python signal_to_noise_pixelConversion.py \
              --selectedRun %s --nCountsPerPhoton %d \
              --ellipse_multiplicative_factor %f --precisionFactor %d \
              --halfWidth %d --label %s'
              %(runNumber, nCountsPerPhoton, 
                ellipse_multiplicative_factor, precisionFactor, 
                halfwidth, label))



label = "_100_lattices"    
flag = 0
if flag == 1:
    os.system('python signal_to_noise_pixelConversion.py \
              --selectedRun %s --nCountsPerPhoton %d \
              --ellipse_multiplicative_factor %f --precisionFactor %d \
              --halfWidth %d --label %s'
              %(runNumber, nCountsPerPhoton, 
                ellipse_multiplicative_factor, precisionFactor, 
                halfwidth, label))


   
flag = 0
if flag == 1:
    os.system('python signal_to_noise_pixelConversion_compare.py --selectedRun %s'%runNumber)



### CC_HALF VS Q, WITH DIFFERENT N LATTICES ###

N_lattices = 10
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d  \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f \
               --precision_factor %d --ellipse_multiplicative_factor %f --N_lattices %d'
               %(runNumber, highResLimit,
                 halfwidth, imagesDirectoryName,
                 geometryFile, referenceCellSize, hmax, kmax,
                 tiltAngle, approximateWavelength, detectorDistance, pixelSize,
                 precisionFactor, ellipse_multiplicative_factor, N_lattices))
                
  
  
N_lattices = 100
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d  \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f \
               --precision_factor %d --ellipse_multiplicative_factor %f --N_lattices %d'
               %(runNumber, highResLimit,
                 halfwidth, imagesDirectoryName,
                 geometryFile, referenceCellSize, hmax, kmax,
                 tiltAngle, approximateWavelength, detectorDistance, pixelSize,
                 precisionFactor, ellipse_multiplicative_factor, N_lattices))
    
    

N_lattices = 586
flag = 0
if flag == 1:
    os.system('python imageSums_displaced_modules_CChalf.py \
               --selectedRun %s --resolutionLimit %f \
               --halfWidth %d --imagesDirectoryName %s  \
               --geometryFile %s --nominalCell %f --hmax %d --kmax %d  \
               --tiltAngle %f --wavelength %f --detectorDistance %f --pixelSize %f \
               --precision_factor %d --ellipse_multiplicative_factor %f --N_lattices %d'
               %(runNumber, highResLimit,
                 halfwidth, imagesDirectoryName,
                 geometryFile, referenceCellSize, hmax, kmax,
                 tiltAngle, approximateWavelength, detectorDistance, pixelSize,
                 precisionFactor, ellipse_multiplicative_factor, N_lattices))


 
flag = 1
if flag == 1:
    os.system('python imageSums_displaced_modules_calculateCC.py --selectedRun %s --cellSize %f'%(runNumber, referenceCellSize))