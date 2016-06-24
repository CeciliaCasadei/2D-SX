# -*- coding: utf-8 -*-
"""
@author: casadei_c
SET selectionFlag ATTRIBUTE TO 1 FOR SELECTED IMAGES.
FOR EACH SELECTED IMAGE, EXTRACT EXPERIMENTAL INFO FROM peaks.txt FILE
(e.g. WAVELENGTH, APPROXIMATE EXPERIMENTAL PEAK POSITIONS AND INTENSITIES....)
AND UPDATE diffractionImage OBJECTS.
"""
import getopt
import sys
import os
import pickle
import h5py
import numpy

import nPeaksDistributionHistogram


def extractExpInfoFunction(myArguments): 
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "selectedImageList=", "peaksFile=", "geometryFile=", "pixelSize="])
    except getopt.GetoptError:
        print 'Usage: python extractExpInfo.py --runNumber <runNumber> --selectedImageList <selectedImageList> --peaksFile <peaksFile> --geometryFile <geometryFile> --pixelSize <pixelSize>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python extractExpInfo.py --runNumber <runNumber> --selectedImageList <selectedImageList> --peaksFile <peaksFile> --geometryFile <geometryFile> --pixelSize <pixelSize>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--selectedImageList":
            selectedImageList = value
        elif option == "--peaksFile":
            peaksFile = value
        elif option == "--geometryFile":
            geometryFile = value
        elif option == "--pixelSize":
            pixelSize = float(value)
            
    # LOAD diffractionImage OBJECTS
    imagesDictionaryFile = './Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber)
    if not os.path.exists(imagesDictionaryFile):
        print 'Data file %s not found.'%imagesDictionaryFile
        return   
    else:
        openPkl = open(imagesDictionaryFile, 'rb')
        myData = pickle.load(openPkl)
        openPkl.close()
        
        # SET IMAGE SELECTION FLAG 
        if not os.path.exists(selectedImageList):
            print 'List file %s not found.'%selectedImageList
            return
        else:
            f = open(selectedImageList, 'r')
            linesList = list(f)
            f.close()
            nSelectedImages = 0
            for i in linesList:
                mySplitLine = i.split()
                myImageName = mySplitLine[1]
                myData['%s'%myImageName].setSelectionFlag(1)
                nSelectedImages = nSelectedImages + 1
            if nSelectedImages == 0:
                print '*** No images were selected! ***'
                return
                
            # EXTRACT GEOMETRY
            if not os.path.exists(geometryFile):
                print 'File %s not found'%geometryFile
                return
            geometryFile = h5py.File(geometryFile, 'r')
            xGeometry = geometryFile['/x']   ### float32 ###
            xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
            yGeometry = geometryFile['/y']   ### float32 ###
            yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
        
    
            # FOR EACH SELECTED IMAGE, EXTRACT EXPERIMENTAL INFO FROM peaks.txt FILE.
            nPeaksList = []                                 
            for i,j in myData.items():
                if  myData['%s'%i].selectionFlag == 1:
                    myData['%s'%i].readPeaksFile(peaksFile, pixelSize, xGeometry_np, yGeometry_np)
                    try:
                        nPeaksList.append(myData['%s'%i].nPeaks)
                    except:
                        continue

            # SAVE UPDATED diffractionImage OBJECTS    
            f = open(imagesDictionaryFile, 'wb') 
            pickle.dump(myData, f)
            f.close()
        
            nPeaksDistributionHistogram.nPeaksDistributionHistogramFunction(runNumber, nPeaksList)
    
    
if __name__ == "__main__":
    print "\n**** CALLING extractExpInfo ****"
    extractExpInfoFunction(sys.argv[1:])