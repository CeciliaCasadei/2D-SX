# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 11:53:56 2016
@author: casadei_c
MEMBER OF diffractionImage CLASS.
FOR EACH DIFFRACTION IMAGE,
EXTRACT EXPERIMENTAL INFORMATION FROM peaks.txt FILE.
GENERATE NEW ATTRIBUTES OF diffractionImage OBJECT:
self.photonEnergy
self.wavelength
self.nPeaks
self.peaksMatrix [x, y, I] -  where x, y are column and row index in unassembled matrix.
self.orderedPeaksMatrix [x, y, I, r, phi, validity] - where x, y are coordinates in pxls, wrt beam center.
"""
import os
import numpy

def readPeaksFileFunction(self, peaksFile, pixelSize, xGeometry_np, yGeometry_np):
    if not os.path.exists(peaksFile):
        print 'File %s not found.'%peaksFile
        return
    else:
        fPeaks = open(peaksFile, 'r')
        expInfoList = list(fPeaks) 
        fPeaks.close()
        
        i = 0
        peak_x_raw_vector = []
        peak_y_raw_vector = []
        I_vector = []
        while i < len(expInfoList):
            listItem = expInfoList[i]
            if '%s'%self.fileName.strip('.h5') in listItem:
                splitLine = listItem.split()
                
                photonEnergy = splitLine[2].strip(',')
                self.photonEnergy = float(photonEnergy)
                
                wavelength = splitLine[3].strip(',')
                self.wavelength = float(wavelength)
                
                peak_x_raw = splitLine[6].strip(',')
                peak_x_raw_vector.append(float(peak_x_raw))
                
                peak_y_raw = splitLine[7].strip(',')
                peak_y_raw_vector.append(float(peak_y_raw))
                
                I = splitLine[12].strip(',')
                I_vector.append(float(I))               
            i = i+1
        self.nPeaks = len(I_vector)  
         
        peaksMatrix = numpy.zeros((self.nPeaks, 3))
        for i in range(0, self.nPeaks):            
            peaksMatrix[i, 0] = peak_x_raw_vector[i]         # Column index in unassembled matrix
            peaksMatrix[i, 1] = peak_y_raw_vector[i]         # Row index in unassembled matrix
            peaksMatrix[i, 2] = I_vector[i]
        self.peaksMatrix = peaksMatrix
        
        # ORDER PEAKS ACCORDING TO DECREASING INTENSITY
        for i in range(0, self.nPeaks):
            for j in range(i+1, self.nPeaks):
                if I_vector[i] < I_vector[j]:
                    x = peak_x_raw_vector[i]
                    y = peak_y_raw_vector[i]
                    Int = I_vector[i]
                    peak_x_raw_vector[i] = peak_x_raw_vector[j]
                    peak_y_raw_vector[i] = peak_y_raw_vector[j]
                    I_vector[i] = I_vector[j]
                    peak_x_raw_vector[j] = x
                    peak_y_raw_vector[j] = y
                    I_vector[j] = Int
            
        # PEAKS COORDINATE CONVERSION
        orderedPeaksMatrix = numpy.zeros((self.nPeaks, 6))
        for i in range(0, self.nPeaks):
            x = xGeometry_np[int(peak_y_raw_vector[i]), int(peak_x_raw_vector[i])]/pixelSize
            y = yGeometry_np[int(peak_y_raw_vector[i]), int(peak_x_raw_vector[i])]/pixelSize
            orderedPeaksMatrix[i,0]= x                                                              ### pxls, wrt beam center ###
            orderedPeaksMatrix[i,1]= y                                                              ### pxls, wrt beam center ###
            orderedPeaksMatrix[i,2]=I_vector[i]
            orderedPeaksMatrix[i,3] = numpy.sqrt(x**2 + y**2)                                       ### radius, pxls ###
            orderedPeaksMatrix[i,4] = numpy.arccos(orderedPeaksMatrix[i,0]/orderedPeaksMatrix[i,3]) ### between 0 and pi ###
            if orderedPeaksMatrix[i,1] < 0:
                orderedPeaksMatrix[i,4] = 2*numpy.pi - orderedPeaksMatrix[i,4]                      ### Azimuth always in [0, 2pi] ###            
            orderedPeaksMatrix[i,5] = 1                                                             ### VALIDITY ###
        self.orderedPeaksMatrix = orderedPeaksMatrix
        
        # LOGGING
        logFolder = './Output_r%s/ObservedPeaksLists'%self.runNumber
        if not os.path.exists(logFolder):   
            os.mkdir(logFolder)  
        f = open('%s/r%s_img%s_ObservedPeaks.txt'%(logFolder, self.runNumber, self.imageNumber), 'w')
        f.write('Writing self.peaksMatrix of diffractionImage object %s\n'%self.fileName)
        f.write('x = column index in unassembled matrix\ny = row index in unassembled matrix\n')
        f.write('           x             y             I')
        for i in range(0, self.nPeaks):        
            f.write('\n%12.2f  %12.2f  %12.2f'%(self.peaksMatrix[i,0], self.peaksMatrix[i,1], self.peaksMatrix[i,2]))
        f.close()
        f = open('%s/r%s_img%s_ObservedPeaks_Sorted.txt'%(logFolder, self.runNumber, self.imageNumber), 'w')
        f.write('Writing self.orderedPeaksMatrix of diffractionImage object %s\n'%self.fileName)
        f.write('x = pxls, wrt beam\ny = pxls, wrt beam\n')
        f.write('           x             y             I             r    azimuth  validity')
        for i in range(0, self.nPeaks):
            f.write('\n%12.2f  %12.2f  %12.2f  %12.2f %10.2f%10d'%(self.orderedPeaksMatrix[i,0], self.orderedPeaksMatrix[i,1], self.orderedPeaksMatrix[i,2], self.orderedPeaksMatrix[i,3], self.orderedPeaksMatrix[i,4], self.orderedPeaksMatrix[i,5]))
        f.close()