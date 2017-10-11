# -*- coding: utf-8 -*-
import getopt
import joblib
import numpy
import sys
import os

def scalingFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["cellSize="])
    except getopt.GetoptError:
        print 'Usage: python model_applyScales_anisotropic.py --cellSize <cellSize>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_applyScales_anisotropic.py --cellSize <cellSize>'
            sys.exit()
        elif option == "--cellSize":
            cellSize = float(value)

    #################################################################################################
    directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
    reciprocalCellRows = 2 * numpy.pi * directCell.I    
    #################################################################################################
        
    runs = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']
               
    # CALCULATE NORMALIZATION FACTOR
    totalScale = 0
    nLattices = 0    
    for run in runs: 
        outputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%run
        LtoModel_vector = joblib.load('%s/r%s-scalesVsModel_anisotropic/r%s-scalesVsModel_anisotropic.jbl'
                                      %(outputFolder, run, run))

        for scaleFactors in LtoModel_vector:
            scale = scaleFactors[0]
            if not numpy.isnan(scale):
                totalScale = totalScale + scale
                nLattices = nLattices + 1
    avgScale = totalScale/nLattices
    print "AVG SCALE ON ALL RUNS: ", avgScale
            
    # CALCULATE NORMALIZED SCALES
    for run in runs:
        outputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%run
        LtoModel_vector = joblib.load('%s/r%s-scalesVsModel_anisotropic/r%s-scalesVsModel_anisotropic.jbl'
                                      %(outputFolder, run, run))
        normalizedScales = []
        for scaleFactors in LtoModel_vector:
            scale = scaleFactors[0]
            if not numpy.isnan(scale):
                normalizedScale = scale/avgScale
            else:
                normalizedScale = numpy.nan
            normalizedScales.append([normalizedScale, scaleFactors[1], scaleFactors[2]])
        joblib.dump(normalizedScales, '%s/r%s-scalesVsModel_anisotropic/r%s-normalizedScalesVsModel_anisotropic.jbl'
                                       %(outputFolder, run, run))   
          
        # APPLY NORMALIZED SCALES
        print '*******APPLY NORMALIZED SCALES TO RUN %s*******'%run          
        myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'
                              %(outputFolder, run, run))
        print len(normalizedScales)
        print len(myList)  
        scaledLatticesList = []
        for lattice in range(0, len(normalizedScales)):
            factors = normalizedScales[lattice]
            scaleFactor  = factors[0]
            damping_para = factors[1]
            damping_perp = factors[2]
            latticeMatrix = myList[lattice]
            if not numpy.isnan(scaleFactor):
                flag = 1
                scaledLattice = []
                for spot in latticeMatrix: # h k qRod I flag i_unassembled j_unassembled
                    h = spot[0]
                    k = spot[1]
                    qRod = spot[2]
                    I = spot[3]
                    
                    reciprocalVector = [h, k]*reciprocalCellRows
                    q_x = reciprocalVector[0,0]         # A^(-1)
                    q_y = reciprocalVector[0,1]         # A^(-1)        
                    q_2D = numpy.sqrt(q_x**2 + q_y**2)  # A^(-1)
            
                    scaledI = scaleFactor * I * numpy.exp(-damping_para*q_2D**2) * numpy.exp(-damping_perp*qRod**2)
                    scaledSpot = [spot[0], spot[1], spot[2], scaledI, flag, spot[5], spot[6]]
                    scaledLattice.append(scaledSpot)
            else:
                flag = 0
                scaledLattice = []
                for spot in latticeMatrix: # h k qRod I flag i_unassembled j_unassembled
                    scaledSpot = [spot[0], spot[1], spot[2], spot[3], flag, spot[5], spot[6]]
                    scaledLattice.append(scaledSpot)
            scaledLattice = numpy.asarray(scaledLattice, dtype = numpy.float32)
            scaledLatticesList.append(scaledLattice)
        
        if not os.path.exists('./Output_runMergingVsModel/spotsMatricesList-Scaled-r%s'%(run)):
            os.mkdir('./Output_runMergingVsModel/spotsMatricesList-Scaled-r%s'%(run))
        joblib.dump(scaledLatticesList, 
                    './Output_runMergingVsModel/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'
                     %(run, run))

    
if __name__ == "__main__":
    print "\n**** CALLING model_applyScales_anisotropic ****"
    scalingFunction(sys.argv[1:])    