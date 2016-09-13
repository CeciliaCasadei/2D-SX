"""
@author: casadei_c
CHECK CHEETAH PREPROCESSING RESULTS
"""
# INSTRUCTIONS: python checkH5content.py --runNumber <runNumber> --label <label>

# -*- coding: utf-8 -*-
import sys
import getopt
import os
import h5py
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot

def checkH5contentFunction(myArguments):
    
    # DEFAULTS
    runNumber = ''
    label = 'good-modified-11'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "label="])
    except getopt.GetoptError:
        print 'Usage: python checkH5content.py --runNumber <runNumber> --label <label>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python checkH5content.py --runNumber <runNumber> --label <label>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--label":
            label = value
    
    # FOLDERS
    if not os.path.exists('./CheckPreprocessing'):
        os.mkdir('./CheckPreprocessing')
    outputFolder = './CheckPreprocessing/r%s-%s'%(runNumber, label)
    if not os.path.exists(outputFolder):
        os.mkdir(outputFolder)
     
    # IMAGE SELECTION 
    myList = sorted(os.listdir('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1'%(runNumber)))
    filenames = myList[0:20]
    
    # READ GEOMETRY    
    geometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5', 'r')
    xGeometry = geometryFile['/x']   ### float32 ###
    xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
    yGeometry = geometryFile['/y']   ### float32 ###
    yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
    
    fOpen = open('%s/log_r%s_%s.log'%(outputFolder, runNumber, label), 'w')
    for filename in filenames:
        
        # LOAD UNASSEMBLED DATA
        print '\nUnassembled data file: %s'%filename
        unassembledDataFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-images/data1/%s'%(runNumber, filename), 'r')
        unassembledData = unassembledDataFile['/data/rawdata0']                       #### int16 #### 
        assembledData   = unassembledDataFile['/data/assembleddata0'] 
        print unassembledData.dtype
        print assembledData.dtype
        unassembledData = numpy.asarray(unassembledData, dtype=numpy.float32)       #### !!!!! ####
        assembledData   = numpy.asarray(assembledData, dtype=numpy.float32)         #### !!!!! ####
        print unassembledData.shape #### 1480, 1552 ####
        print assembledData.shape   #### 1750, 1750 #### 
    
        # PLOT UNASSEMBLED IMAGE
        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
        matplotlib.pyplot.title('%s'%filename, y=1.05)
        myAxesImageObject = matplotlib.pyplot.imshow(unassembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        matplotlib.pyplot.savefig('%s/unassembledImage_%s.png'%(outputFolder, filename.strip('.h5')))
        matplotlib.pyplot.close()   
                                                             
        # PLOT ASSEMBLED IMAGE       
        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
        matplotlib.pyplot.title('%s'%filename, y=1.05)
        myAxesImageObject = matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        matplotlib.pyplot.savefig('%s/assembledImage_%s.png'%(outputFolder, filename.strip('.h5')))
        matplotlib.pyplot.close()   
        
        # EXTRACT PEAKS COORDINATES FROM CHEETAH peaks.txt FILE   
        newPeaksFile = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r%s-%s/peaks.txt'%(runNumber, label)
        fPeaks = open(newPeaksFile, 'r')
        expInfoList = list(fPeaks) 
        fPeaks.close()
        peak_x_raw = []
        peak_y_raw = []
        lineN = 0
        for myLine in expInfoList:
            if lineN > 0:
                splitLine = myLine.split()
                if splitLine[1].strip(',') == filename.strip('.h5'):
                    x_raw = float(splitLine[6].strip(','))
                    y_raw = float(splitLine[7].strip(','))
                    peak_x_raw.append(x_raw)
                    peak_y_raw.append(y_raw)
                    print '\nSpot area in newUnassembledData:'
                    print unassembledData[int(y_raw)-5:int(y_raw)+5, int(x_raw)-5:int(x_raw)+5]
            lineN = lineN + 1
        peak_x_raw = numpy.asarray(peak_x_raw, dtype=numpy.float32)  # column index in unassembled data matrix
        peak_y_raw = numpy.asarray(peak_y_raw, dtype=numpy.float32)  # row index in unassembled data matrix
        peak_x_raw_converted = []
        peak_y_raw_converted = []
        for i in range(0, peak_x_raw.shape[0]):
            x = xGeometry_np[int(peak_y_raw[i]), int(peak_x_raw[i])]/0.000110 # pxls, wrt beam center
            y = yGeometry_np[int(peak_y_raw[i]), int(peak_x_raw[i])]/0.000110 # pxls, wrt beam center
            peak_x_raw_converted.append(x)
            peak_y_raw_converted.append(y)
    
        # N OF PEAKS FOUND
        print 'Cheetah peaks: %d'%(len(peak_x_raw_converted))
        fOpen.write('%s: n Cheetah peaks: %d\n'%(filename, len(peak_x_raw_converted)))
            
        # PLOT ASSEMBLED DATA AND PEAKS FOUND IN NEW peaks.txt
        peak_x_raw_converted = numpy.asarray(peak_x_raw_converted)
        peak_x_raw_converted = peak_x_raw_converted + 875   ### or 874 ###
        peak_y_raw_converted = numpy.asarray(peak_y_raw_converted)
        peak_y_raw_converted = peak_y_raw_converted + 875   ### or 874 ###
        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
        matplotlib.pyplot.title('%s'%filename, y=1.05)
        myAxesImageObject = matplotlib.pyplot.imshow(assembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        matplotlib.pyplot.scatter(peak_x_raw_converted, peak_y_raw_converted, marker="o", color='w', facecolors='none', s=40)
        matplotlib.pyplot.savefig('%s/cheetahPeaks_image_%s.png'%(outputFolder, filename.strip('.h5')))
        matplotlib.pyplot.close()  
        
        # PLOT UNASSEMBLED DATA AND PEAKS FOUND IN NEW peaks.txt
        myFigureObject  = matplotlib.pyplot.figure(figsize=(40,40), dpi=4*96, facecolor='w',frameon=True)
        matplotlib.pyplot.title('%s'%filename, y=1.05)
        myAxesImageObject = matplotlib.pyplot.imshow(unassembledData, origin='lower', interpolation='nearest', vmin = 0, vmax = 100)
        myFigureObject.colorbar(myAxesImageObject, pad=0.01, fraction=0.0471, shrink=1.00, aspect=20)
        matplotlib.pyplot.scatter(peak_x_raw, peak_y_raw, marker="o", color='w', facecolors='none', s=40)
        matplotlib.pyplot.savefig('%s/cheetahPeaks_unassembled_image_%s.png'%(outputFolder, filename.strip('.h5')))
        matplotlib.pyplot.close()  
    fOpen.close()
    
if __name__ == "__main__":
    print "\n**** CALLING checkH5content ****"
    checkH5contentFunction(sys.argv[1:])