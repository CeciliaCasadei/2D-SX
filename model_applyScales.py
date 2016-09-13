# -*- coding: utf-8 -*-
import sys
import getopt
import joblib
import numpy
import os


def model_applyScalesFunction(myArguments):
    runNumbers = ['0127', '0195', '0196', '0197', '0198', '0199', '0200', '0201']
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder="])
    except getopt.GetoptError:
        print 'Usage: python model_applyScales.py'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyScales.py'
            sys.exit()
    
    # CALCULATE NORMALIZATION FACTOR
    totalScale = 0
    nLattices = 0    
    for runNumber in runNumbers:
        scales = joblib.load('./Output_runMergingVsModel/transformAndScaleToModel_r%s/r%s-scales/r%s-scales.jbl'%(runNumber, runNumber, runNumber))
        print runNumber
        print len(scales)
        for scale in scales:
            if not numpy.isnan(scale):
                totalScale = totalScale + scale
                nLattices = nLattices + 1
    avgScale = totalScale/nLattices
    
    # CALCULATE NORMALIZED SCALES
    for runNumber in runNumbers:
        scales = joblib.load('./Output_runMergingVsModel/transformAndScaleToModel_r%s/r%s-scales/r%s-scales.jbl'%(runNumber, runNumber, runNumber))
        normalizedScales = []
        print runNumber
        print len(scales)
        for scale in scales:
            if not numpy.isnan(scale):
                normalizedScale = scale/avgScale
            else:
                normalizedScale = numpy.nan
            normalizedScales.append(normalizedScale)
            
        joblib.dump(normalizedScales, './Output_runMergingVsModel/transformAndScaleToModel_r%s/r%s-scales/r%s-normalizedScales.jbl'%(runNumber, runNumber, runNumber))
      
    # APPLY NORMALIZED SCALES 
    for runNumber in runNumbers:
        normalizedScales = joblib.load('./Output_runMergingVsModel/transformAndScaleToModel_r%s/r%s-scales/r%s-normalizedScales.jbl'%(runNumber, runNumber, runNumber))
        latticesList = joblib.load('./Output_runMergingVsModel/transformAndScaleToModel_r%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(runNumber, runNumber, runNumber))
        print '****************'        
        print len(normalizedScales)
        print len(latticesList)  
        scaledLatticesList = []
        for lattice in range(0, len(normalizedScales)):
            normalizedScale = normalizedScales[lattice]
            latticeMatrix = latticesList[lattice]
            if not numpy.isnan(normalizedScale):
                flag = 1
                scaledLattice = []
                for spot in latticeMatrix: # h k qRod I flag i_unassembled j_unassembled
                    scaledSpot = [spot[0], spot[1], spot[2], normalizedScale*spot[3], flag, spot[5], spot[6]]
                    scaledLattice.append(scaledSpot)
            else:
                flag = 0
                scaledLattice = []
                for spot in latticeMatrix: # h k qRod I flag i_unassembled j_unassembled
                    scaledSpot = [spot[0], spot[1], spot[2], spot[3], flag, spot[5], spot[6]]
                    scaledLattice.append(scaledSpot)
            scaledLattice = numpy.asarray(scaledLattice, dtype = numpy.float32)
            scaledLatticesList.append(scaledLattice)
        
        if not os.path.exists('./Output_runMergingVsModel/spotsMatricesList-Scaled-r%s'%runNumber):
            os.mkdir('./Output_runMergingVsModel/spotsMatricesList-Scaled-r%s'%runNumber)
        joblib.dump(scaledLatticesList, './Output_runMergingVsModel/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(runNumber, runNumber))
        
     
if __name__ == "__main__":
    print "\n**** CALLING model_applyScales ****"
    model_applyScalesFunction(sys.argv[1:])    