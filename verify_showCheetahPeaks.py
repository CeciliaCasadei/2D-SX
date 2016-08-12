# -*- coding: utf-8 -*-
import getopt
import sys
import os
import pickle
import numpy
import h5py
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot
import warnings

def verify_showCheetahPeaksFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "imageNumber=", "imageFolderName=", "geometryFile="])
    except getopt.GetoptError:
        print 'Usage: python verify_showCheetahPeaks.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName> --geometryFile <geometryFile>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python verify_showCheetahPeaks.py --runNumber <runNumber> --imageNumber <imageNumber> --imageFolderName <imageFolderName> --geometryFile <geometryFile>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
        elif option == "--imageNumber":
            imageNumber = value
        elif option == "--imageFolderName":
            imageFolderName = value
        elif option == "--geometryFile":
            geometryFile = value
    
    # READ GEOMETRY
    if not os.path.exists(geometryFile):
        print 'File %s not found'%geometryFile
        return
    geometryFile = h5py.File(geometryFile, 'r')
    xGeometry = geometryFile['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryFile['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
        
    # DISPLAY ASSEMBLED IMAGE AND CHEETAH PEAKS
    pklFile = './Output_r%s/ExtractExperimentalInfo/r%s_imagesDictionary.pkl'%(runNumber, runNumber)
    if not os.path.exists(pklFile):
        print 'File %s not found'%pklFile
    else:
        openPkl = open(pklFile, 'rb')
        imageObjects = pickle.load(openPkl)
        openPkl.close()
        if not os.path.exists(imageFolderName):
            print 'File %s not found'%imageFolderName
        else:
            for imageKey, imageObject in imageObjects.items():
                if imageObject.runNumber == runNumber and imageObject.imageNumber == imageNumber:
                    try:
                        peaksMatrix = imageObject.peaksMatrix
                        x_vector = []
                        y_vector = []
                        x_converted_vector = []
                        y_converted_vector = []
                        for peak in peaksMatrix:
                            x = peak[0] # column index in unassembled file
                            y = peak[1] # row index in unassembled file
                            x_vector.append(x)
                            y_vector.append(y)
                            x_converted = xGeometry_np[int(y), int(x)]/0.000110 # x coo. in pxls, wrt beam 
                            y_converted = yGeometry_np[int(y), int(x)]/0.000110 # y coo. in pxls, wrt beam 
                            x_converted = x_converted + 875
                            y_converted = y_converted + 875
                            x_converted_vector.append(x_converted)
                            y_converted_vector.append(y_converted)
                        unassembledDataFile = '%s/%s'%(imageFolderName, imageObject.fileName)
                        unassembledData = h5py.File(unassembledDataFile, 'r')
                        
                        assembledData0 = unassembledData['/data/assembleddata0'] #### (1750, 1750) ####
                        warnings.filterwarnings("ignore")
                        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
                        matplotlib.pyplot.title('%s\nRun: %s - Number: %s'%(imageObject.fileName, imageObject.runNumber, imageObject.imageNumber), y=1.05)
                        myAxesImageObject = matplotlib.pyplot.imshow(assembledData0, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
                        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                        matplotlib.pyplot.scatter(x_converted_vector, y_converted_vector, marker="o", color='w', facecolors='none', s=40)
                        if not os.path.exists('./Output_r%s/CheetahPeaks'%runNumber):
                            os.mkdir('./Output_r%s/CheetahPeaks'%runNumber)
                        matplotlib.pyplot.savefig('./Output_r%s/CheetahPeaks/run%s_image%s_assembled.png'%(runNumber, imageObject.runNumber, imageObject.imageNumber))
                        matplotlib.pyplot.close()  
                        
                        rawData0 = unassembledData['/data/rawdata0'] #### () ####
                        warnings.filterwarnings("ignore")
                        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
                        matplotlib.pyplot.title('%s\nRun: %s - Number: %s'%(imageObject.fileName, imageObject.runNumber, imageObject.imageNumber), y=1.05)
                        myAxesImageObject = matplotlib.pyplot.imshow(rawData0, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
                        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
                        matplotlib.pyplot.scatter(x_vector, y_vector, marker="o", color='w', facecolors='none', s=40)        
                        matplotlib.pyplot.savefig('./Output_r%s/CheetahPeaks/run%s_image%s_unassembled.png'%(runNumber, imageObject.runNumber, imageObject.imageNumber))
                        matplotlib.pyplot.close()   
                        
                    except:
                        print 'Exception'
                        return
                    
                    
                    
if __name__ == "__main__":
    print "\n**** SHOW CHEETAH PEAKS ****"
    verify_showCheetahPeaksFunction(sys.argv[1:])