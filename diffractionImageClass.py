# -*- coding: utf-8 -*-
class diffractionImage: 
    
    def __init__(self, fileName, runNumber, imageNumber, tiltAngle, selectionFlag):
        self.fileName = fileName
        self.runNumber = runNumber
        self.imageNumber = imageNumber
        self.tiltAngle = tiltAngle
        self.selectionFlag = selectionFlag
        # diffractionImage objects generated in storeImageObjects.py
    

        
    def displayImage(self, folderName):
        import displayImage
        displayImage.displayImageFunction(self, folderName)
        # Display the assembled image associated to each object.
    

     
    def setSelectionFlag(self, selectionFlag):
        self.selectionFlag = selectionFlag
        # Set image selection flag in extractExpInfo.py
    

    
    def readPeaksFile(self, peaksFile, pixelSize, xGeometry_np, yGeometry_np):
        import readPeaksFile
        readPeaksFile.readPeaksFileFunction(self, peaksFile, pixelSize, xGeometry_np, yGeometry_np)
        # Extract experimental info from peaks.txt file.
        # Used in extractExpInfo.py
        # Generate new attributes of diffractionImage object:
        # self.nPeaks
        # self.photonEnergy
        # self.wavelength
        # self.peaksMatrix [x, y, I] -  where x, y are column and row index in unassembled matrix.
        # Sort experimental peaks in order of decreasing intensity
        # Convert peak coordinates from unassembled matrix indices
        # to x, y in pxls, wrt beam center and radius, azimuth
        # Store result in self.orderedPeaksMatrix 
        # self.orderedPeaksMatrix [x, y, I, r, phi, validity] - where x, y are coordinates in pxls, wrt beam center


    
    def indexingFunction(self, detectorDistance, pixelSize, 
                         radialTolerance, pixelTolerance, azimuthTolerance,
                         minNofPeaksPerLattice, maxNofPeaksPerImage, 
                         referenceCellSize,
                         geometryFile):
        if self.nPeaks <= maxNofPeaksPerImage:
            import indexing
            indexing.indexingFunction(self, detectorDistance, pixelSize,
                                      radialTolerance, pixelTolerance, azimuthTolerance,
                                      minNofPeaksPerLattice, 
                                      referenceCellSize,
                                      geometryFile)
        else:
            print 'Image %s - %s:\tnumber of peaks above threshold.'%(self.imageNumber.zfill(5), self.fileName)
        # Identify and index multiple lattices in one image.
        # Generate objects of the class Lattice.


        
    def getPredictedPattern(self, detectorDistance, pixelSize, cellSize, trialInPlaneAngles):
        import getPredictedPattern
        getPredictedPattern.getPredictedPatternFunction(self, detectorDistance, pixelSize, cellSize, trialInPlaneAngles)
        # Calculate trial predicted patterns in 256 in-plane lattice orientations
        # Generate new attribute:
        # self.referencePredictedPattern: dictionary containing 256 predicted patterns.
        #                                 [h k qx qy dMin q azimuth rotatedAzimuth detectorAzimuth diffractionAngle detectorRadius qRod LPfactor]
     

                                                   
    def plotIndexedExperimentalPeaks(self, detectorDistance, 
                                     pixelSize, resolutionRadii):
        import plotIndexedExperimentalPeaks
        plotIndexedExperimentalPeaks.plotIndexedExperimentalPeaksFunction(self, detectorDistance, 
                                                                          pixelSize, resolutionRadii) 
        # Plot indexing results
        # Plot experimental peaks (from peaks.txt) and indexed calculated pattern(s)
    

                                                                      
    def plotRefinedLattices(self, resolutionRadii):
        import plotRefinedIndexedImages
        plotRefinedIndexedImages.plotRefinedIndexedImagesFunction(self, resolutionRadii)
        # Plot refinement results
        # Plot experimental peaks (from peaks.txt) and indexed calculated pattern(s)
        # after lattice orientation and cell size refinement.
        
    
    def plotRefinedLattices_imageOverlap(self, resolutionRadii):
        import plotRefinedIndexedImages
        plotRefinedIndexedImages.plotRefinedIndexedImagesFunction_imageOverlap(self, resolutionRadii)