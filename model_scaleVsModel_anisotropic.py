# -*- coding: utf-8 -*-
import getopt
import joblib
import time
import numpy
import sys
import os
import scipy.optimize

import scaling_matchToModel_anisotropic

def Correlate(x1, x2):
    x1Avg = numpy.average(x1)
    x2Avg = numpy.average(x2)
    numTerm = numpy.multiply(x1-x1Avg, x2-x2Avg)
    num = numTerm.sum()
    resX1Sq = numpy.multiply(x1-x1Avg, x1-x1Avg)
    resX2Sq = numpy.multiply(x2-x2Avg, x2-x2Avg)
    den = numpy.sqrt(numpy.multiply(resX1Sq.sum(), resX2Sq.sum()))
    CC = num/den
    return CC
    
def toBeMinimized(x, *p):
    scaleFactor = x[0]
    TFactor_para = x[1]
    TFactor_perp = x[2]
    
    I_model, I, q_2D, q_rod = p
    I_model = numpy.array(I_model)
    I = numpy.array(I)
    q_2D = numpy.array(q_2D)
    q_rod = numpy.asarray(q_rod)
    
    scaled_I_vector = (scaleFactor * I * 
                       numpy.exp(-TFactor_para*q_2D**2) * 
                       numpy.exp(-TFactor_perp*q_rod**2))
    intensityError = (I_model - scaled_I_vector)**2
    intensityError = intensityError.sum()
    return intensityError

def scalingFunction(myArguments):
    
    str_1 = '--runNumber <runNumber> --CC_threshold <CC_threshold>'
    str_2 = '--cellSize <cellSize> --resolution_3D <resolution_3D>'

    # DEFAULTS
    deltaQrodThreshold = 0.001
    n_minThreshold = 100
    CC_threshold = 0.94
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["runNumber=", "CC_threshold=", 
                                               "cellSize=",  "resolution_3D="])
    except getopt.GetoptError:
        print 'Usage: python model_scaleVsModel_anisotropic.py %s'%(str_1, str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_scaleVsModel_anisotropic.py %s'%(str_1, str_2)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value
        elif option == "--CC_threshold":
            CC_threshold = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
        elif option == "--resolution_3D":
            resolution_3D = float(value)
            print resolution_3D
    
    inputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%runNumber
    
    # LOAD MODEL: h k qRod I        
    lattice_model = joblib.load('./Output_runMerging/Shannon_sampling/model/lattice_model.jbl')
    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
           
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag i_unassembled j_unassembled
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'
                          %(inputFolder, runNumber, runNumber))
    nLattices = len(myList)
    print 'Run: %s Total n of lattices: %d'%(runNumber, nLattices)
    
    # SCALE LATTICES WITH RESPECT TO MODEL    
    LtoModel_vector = []
    nScaled = 0
    
    startTime = time.time()  
    fOpen = open('%s/r%s_scalingVsModel_anisotropic.txt'%(inputFolder, runNumber), 'w')
    fOpen.write('Total: %s lattices\n'%nLattices)
    fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
    fOpen.write('Min pairs n: %d\n'%n_minThreshold)
    fOpen.write('CC threshold: %f\n'%CC_threshold)
    
    for firstNeighbor in range(0, nLattices):
        spots1stN = myList[firstNeighbor]  
        if spots1stN[0, 4] == 0:
            refinedScale  = numpy.nan
            refinedT_para = numpy.nan
            refinedT_perp = numpy.nan
            print 'Lattice %s: Non oriented\n'%firstNeighbor
            fOpen.write('Lattice %s: Non oriented\n'%firstNeighbor)            
        else:        
            print 'Lattice %s'%firstNeighbor                

            n_pairs, \
            I1, \
            I2, \
            qs_2D, \
            qs_rod = scaling_matchToModel_anisotropic.matching(lattice_model, 
                                                               spots1stN, 
                                                               deltaQrodThreshold,
                                                               cellSize,
                                                               resolution_3D)
            if n_pairs < n_minThreshold:
                refinedScale  = numpy.nan
                refinedT_para = numpy.nan
                refinedT_perp = numpy.nan
                fOpen.write('Lattice %s: n_pairs model-lattice below threshold.\n'
                             %firstNeighbor)
            else:  
                # I1 = model
                # I2 = lattice
                ##########################################################################
                I1 = numpy.asarray(I1)
                I2 = numpy.asarray(I2)
                qs_2D = numpy.asarray(qs_2D)
                qs_rod = numpy.asarray(qs_rod)
                I1 = I1.flatten()
                I2 = I2.flatten()
                qs_2D = qs_2D.flatten()
                qs_rod = qs_rod.flatten()
                print I1.shape, I2.shape, qs_2D.shape, qs_rod.shape
                
                X_0 = 1
                B_para_0 = 0
                B_perp_0 = 0
                xStart = [X_0, B_para_0, B_perp_0]
                fixedParas = (I1, I2, qs_2D, qs_rod)
                result = scipy.optimize.minimize(toBeMinimized, 
                                                 xStart, 
                                                 args = fixedParas, 
                                                 method = "Powell", 
                                                 options = {'xtol': 0.0001, 'ftol': 0.0001})
                print result.message
                ### RESULTS ###
                refinedScale = result.x[0]
                refinedT_para = result.x[1]
                refinedT_perp = result.x[2]
                print refinedScale, refinedT_para, refinedT_perp
                
                scaled_I2 = (refinedScale * I2 * 
                             numpy.exp(-refinedT_para*qs_2D**2) * 
                             numpy.exp(-refinedT_perp*qs_rod**2))
                CC = Correlate(I1, scaled_I2)
                print 'N points: %d'%len(I1)
                print 'CC: %.4f'%CC
                ##########################################################################
                if CC >= CC_threshold:
                    nScaled = nScaled + 1

                    fOpen.write('Lattice %s, Scale (lattice to model): %.3f, T (lattice to model): %.3f %.3f\n'
                                %(firstNeighbor, refinedScale, refinedT_para, refinedT_perp))
                    print ('Lattice %s, Scale (lattice to model): %.3f, T (lattice to model): %.3f %.3f\n'
                            %(firstNeighbor, refinedScale, refinedT_para, refinedT_perp))
                else:
                    refinedScale  = numpy.nan
                    refinedT_para = numpy.nan
                    refinedT_perp = numpy.nan
                    fOpen.write('Lattice %s, CC below threshold.\n'%firstNeighbor) 
                    print 'Lattice %s, CC below threshold.\n'%firstNeighbor
               
        LtoModel_vector.append([refinedScale, refinedT_para, refinedT_perp]) # lattice to model
    scaledFraction = float(nScaled)/nLattices
    fOpen.write('\nFraction of scaled lattices: %.2f.'%scaledFraction)    
    runTime = time.time() - startTime
    fOpen.write('\nIt took: %.1f s'%runTime)
    fOpen.close
    print 'Total n of lattices: %d'%len(LtoModel_vector)  
    print 'Fraction of scaled lattices: %.2f.'%scaledFraction
 
    if not os.path.exists('%s/r%s-scalesVsModel_anisotropic'%(inputFolder, runNumber)):
        os.mkdir('%s/r%s-scalesVsModel_anisotropic'%(inputFolder, runNumber))           
    joblib.dump(LtoModel_vector, '%s/r%s-scalesVsModel_anisotropic/r%s-scalesVsModel_anisotropic.jbl'
                                  %(inputFolder, runNumber, runNumber))
           
    
if __name__ == "__main__":
    print "\n**** CALLING model_scaleVsModel_anisotropic ****"
    scalingFunction(sys.argv[1:])    