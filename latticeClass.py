# -*- coding: utf-8 -*-
class Lattice:
    
    def __init__(self, fileName, imageNumber, runNumber, tiltAngle,
                 inPlaneRotation, wavelength, pixelSize, detectorDistance, 
                 nMatchedPeaks, referenceCellSize, indexedPeaksTable,
                 latticeNumber):
        self.fileName = fileName
        self.imageNumber = imageNumber
        self.runNumber = runNumber
        self.tiltAngle = tiltAngle
        self.inPlaneRotation = inPlaneRotation
        self.wavelength = wavelength
        self.pixelSize = pixelSize
        self.detectorDistance = detectorDistance
        self.nMatchedPeaks = nMatchedPeaks
        self.referenceCellSize = referenceCellSize
        self.indexedPeaksTable = indexedPeaksTable
        self.latticeNumberInImage = latticeNumber
        # Lattice objects are generated in indexing.py (member of diffractionImage class) 
        # self.indexedPeaksTable:
        # [h k 
        #  experimentalRadius experimentalAzimuth experimentalI 
        #  radialDelta azimuthDelta 
        #  experimentalPeakN predictedPeakN 
        #  predictedRadius predictedAzimuth]
        
        
        
    def refineCellSizeAndOrientation(self, 
                                     sizeRefinementSteps, orientationRefinementSteps, 
                                     folderName, RLfolderName):
        import refineCellSizeAndOrientation
        refinedParameters = refineCellSizeAndOrientation.refineCellSizeAndOrientation(self, 
                                                                                      sizeRefinementSteps, orientationRefinementSteps, 
                                                                                      folderName, RLfolderName)
        return refinedParameters
        # Called by orientationAndCellRefinement.py 
        # Given a list of cell size refinement steps and a list of in-plane orientation refinement steps,
        # the lattice error is calculated in each grid point of the 2D refinement parameters space and minimized.
          
          
          
    def setRefinedCellSize(self, refinedCellSize):
        self.refinedCellSize = refinedCellSize
    def setRefinedInPlaneOrientation(self, refinedOrientation):
        self.refinedInPlaneOrientation = refinedOrientation
    def setRefinedPattern(self, refinedPattern):
        self.refinedPredictedPattern = refinedPattern
    def setLatticeErrorMatrix(self, latticeErrorMatrix):
        self.latticeErrorMatrix = latticeErrorMatrix
    def setLatticeError(self, avgMinError):
        self.avgLatticeError = avgMinError
        ### Refined Lattice attributes set in orientationAndCellRefinement.py ###
         
         
         
    def logRefinedPatternInfo(self, folderName):
        import logRefinedPatternInfo
        logRefinedPatternInfo.logRefinedPatternInfoFunction(self, folderName)
        # Called by orientationAndCellRefinement.py
        # Output a text file reporting the refined predicted pattern data.



    def plotLatticeErrorMatrix(self, filename):
        import plotLatticeErrorMatrix
        plotLatticeErrorMatrix.plotLatticeErrorMatrixFunction(self, filename)
        # Called by orientationAndCellRefinement.py
        # Plot the attribute self.LatticeErrorMatrix
    


    def refinementBehaviourPlot(self, folderName, myString):
        import refinementBehaviourPlot
        refinementBehaviourPlot.refinementBehaviourPlotFunction(self, folderName, myString)
        # Called by processingAndIntegration
        # Plot behaviour of in-plane orientation, cell size, center coordinates,
        # n of detected spots and lattice error
        # during refinement iterations in image processing