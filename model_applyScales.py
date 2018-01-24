# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os


def model_applyScalesFunction(myArguments):
    runNumbers = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h")
    except getopt.GetoptError:
        print 'Usage: python model_applyScales.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyScales.py'
            sys.exit()
    
    baseFolder = './Output_runMergingVsModel'
    # CALCULATE NORMALIZATION FACTOR
    totalScale = 0
    nLattices = 0    
    for runNumber in runNumbers:
        folder = '%s/transformAndScaleToModel_r%s'%(baseFolder, runNumber)
        scales = joblib.load('%s/r%s-scales/r%s-scales.jbl'%(folder, 
                                                             runNumber, 
                                                             runNumber))
        print 'RUN ', runNumber
        print len(scales), ' LATTICES'
        for scale in scales:
            if not numpy.isnan(scale):
                totalScale = totalScale + scale
                nLattices = nLattices + 1
    avgScale = totalScale/nLattices
    
    # CALCULATE NORMALIZED SCALES
    for runNumber in runNumbers:
        folder = '%s/transformAndScaleToModel_r%s'%(baseFolder, runNumber)
        scales = joblib.load('%s/r%s-scales/r%s-scales.jbl'%(folder, 
                                                             runNumber, 
                                                             runNumber))
        normalizedScales = []
        print 'RUN ', runNumber
        print len(scales), ' LATTICES'
        for scale in scales:
            if not numpy.isnan(scale):
                normalizedScale = scale/avgScale
            else:
                normalizedScale = numpy.nan
            normalizedScales.append(normalizedScale)
            
        joblib.dump(normalizedScales, 
                    '%s/r%s-scales/r%s-normalizedScales.jbl'%(folder, 
                                                              runNumber, 
                                                              runNumber))
      
    # APPLY NORMALIZED SCALES 
    for runNumber in runNumbers:
        folder = '%s/transformAndScaleToModel_r%s'%(baseFolder, runNumber)
        normalizedScales = joblib.load('%s/r%s-scales/r%s-normalizedScales.jbl'
                                        %(folder, runNumber, runNumber))
        folder_t = '%s/spotsMatricesList-Transformed-r%s'%(folder, runNumber)
        latticesList = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                                    %(folder_t, runNumber))
        print '****************'        
        print len(normalizedScales), ' LATTICES' 
        print len(latticesList), ' LATTICES'  
        scaledLatticesList = []
        for lattice in range(0, len(normalizedScales)):
            normalizedScale = normalizedScales[lattice]
            latticeMatrix = latticesList[lattice]
            if not numpy.isnan(normalizedScale):
                flag = 1
                scaledLattice = []
                for spot in latticeMatrix: 
                    scaledSpot = [spot[0],                     # h
                                  spot[1],                     # k
                                  spot[2],                     # qRod
                                  normalizedScale*spot[3],     # I
                                  flag,                        # flag
                                  spot[5],                     # i
                                  spot[6]]                     # j
                    scaledLattice.append(scaledSpot)
            else:
                flag = 0
                scaledLattice = []
                for spot in latticeMatrix: 
                    scaledSpot = [spot[0], 
                                  spot[1], 
                                  spot[2], 
                                  spot[3], 
                                  flag, 
                                  spot[5], 
                                  spot[6]]
                    scaledLattice.append(scaledSpot)
            scaledLattice = numpy.asarray(scaledLattice, dtype = numpy.float32)
            scaledLatticesList.append(scaledLattice)
        
        folder_s = '%s/spotsMatricesList-Scaled-r%s'%(baseFolder, runNumber)
        if not os.path.exists(folder_s):
            os.mkdir(folder_s)
        joblib.dump(scaledLatticesList,
                    '%s/r%s_scaledSpotsMatricesList.jbl'%(folder_s, runNumber))
        
     
if __name__ == "__main__":
    print "\n**** CALLING model_applyScales ****"
    model_applyScalesFunction(sys.argv[1:])    