# -*- coding: utf-8 -*-
# diffractionSpot objects are generated in unassembledImageProcessing
class diffractionSpot:
    
    def __init__(self, n, h, k):
        self.n = n
        self.h = h
        self.k = k


        
    def setBoxMatrix(self, boxMatrix):
        self.boxMatrix = boxMatrix
        
    def setMaskBoxMatrix(self, maskBoxMatrix):
        self.maskBoxMatrix = maskBoxMatrix
        
    def setBgSubtractedBoxMatrix(self, bgSubtractedBoxMatrix):
        self.bgSubtractedBoxMatrix = bgSubtractedBoxMatrix
        
    def setPeakDetectionMask(self, peakDetectionMask):
        self.peakDetectionMask = peakDetectionMask
        
    def setBoxLimits(self, xLeft, yDown):
        self.xLeft = xLeft
        self.yDown = yDown        

    def setCoM(self, xCoM, yCoM): 
        self.xCoM = xCoM
        self.yCoM = yCoM
     
    def setConnectedI(self, connectedI):
        self.connectedI = connectedI
        
    def setDistanceFromPrediction(self, distanceFromPrediction):
        self.distanceFromPrediction = distanceFromPrediction
        
        
        
    def setBoxIndices(self, i, j):
        self.i = i
        self.j = j
        
        
        
    def setFinalBoxIndices(self, i, j):
        self.iFinal = i
        self.jFinal = j
        
        
        
    def integrateSpot(self):
        import integrateSpot
        integratedIntensity = integrateSpot.integrate(self)
        self.integratedIntensity = integratedIntensity
        
        
        
    def spotProcessingPlot(self, runNumber, imageNumber, latticeNumberInImage, boxWidth):
        import processingPlot
        processingPlot.processingPlotFunction(self, runNumber, imageNumber, latticeNumberInImage, boxWidth)
        