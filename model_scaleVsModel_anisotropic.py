# -*- coding: utf-8 -*-
import getopt
import joblib
import time
import numpy
import sys
import os
import scipy.optimize
import matplotlib.pyplot

import scaling_matchToModel_anisotropic
import correlate
    
def toBeMinimized_square(x, *p):
    scaleFactor  = x[0]
    TFactor_para = x[1]
    TFactor_perp = x[2]
    
    I_model, I, q_2D, q_rod = p
    I_model = numpy.asarray(I_model)
    I = numpy.asarray(I)
    q_2D = numpy.asarray(q_2D)
    q_rod = numpy.asarray(q_rod)
    
    scaled_I_vector = (scaleFactor * I * 
                       numpy.exp(-TFactor_para*q_2D**2) * 
                       numpy.exp(-TFactor_perp*q_rod**2))
    intensityError = (I_model - scaled_I_vector)**2
    intensityError = intensityError.sum()
    return intensityError

def toBeMinimized_linear(x, *p):
    scaleFactor  = x[0]
    TFactor_para = x[1]
    TFactor_perp = x[2]
    
    I_model, I, q_2D, q_rod = p
    I_model = numpy.asarray(I_model)
    I = numpy.asarray(I)
    q_2D = numpy.asarray(q_2D)
    q_rod = numpy.asarray(q_rod)
    
    scaled_I_vector = (scaleFactor * I * 
                       numpy.exp(-TFactor_para*q_2D**2) * 
                       numpy.exp(-TFactor_perp*q_rod**2))
    intensityError = abs(I_model - scaled_I_vector)
    intensityError = intensityError.sum()
    return intensityError
    
def scalingFunction(myArguments):
    
    str_1 = '--runNumber <runNumber> --CC_threshold <CC_threshold>'
    str_2 = '--cellSize <cellSize> --resolution_3D <resolution_3D>'

    # DEFAULTS
    deltaQrodThreshold = 0.001/2
    n_minThreshold = 30 #100
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                              ["runNumber=", "CC_threshold=", 
                                               "cellSize=",  "resolution_3D="])
    except getopt.GetoptError:
        print 'Usage: python model_scaleVsModel_anisotropic.py %s'%(str_1, 
                                                                    str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_scaleVsModel_anisotropic.py %s'%(str_1, 
                                                                        str_2)
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
    
    inputPath = './Output_runMergingVsModel'
    inputFolder = '%s/transformAndScaleToModel_r%s'%(inputPath, runNumber)
    figuresFolder = '%s/Figures'%inputFolder
    if not os.path.exists(figuresFolder):
        os.mkdir(figuresFolder)
    
    # LOAD MODEL: h k qRod I 
    lattice_model_path = './Output_runMerging/Shannon_sampling/model'       
    lattice_model = joblib.load('%s/lattice_model.jbl'%lattice_model_path)
    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
           
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag i j
    listPath = '%s/spotsMatricesList-Transformed-r%s'%(inputFolder, runNumber)
    myList = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                          %(listPath, runNumber))
    nLattices = len(myList)
    print 'Run: %s Total n of lattices: %d'%(runNumber, nLattices)
    
    # SCALE LATTICES WITH RESPECT TO MODEL    
    LtoModel_vector = []
    nScaled = 0
    
    startTime = time.time()  
    fOpen = open('%s/r%s_scalingVsModel_anisotropic.txt'%(inputFolder, 
                                                          runNumber), 
                                                          'w')
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
            print 'L %s: Non oriented\n'%firstNeighbor
            fOpen.write('L %s: Non oriented\n'%firstNeighbor)            
        else:        
            print 'L %s'%firstNeighbor                

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
                fOpen.write('L %s: n_pairs model-lattice below threshold.\n'
                             %firstNeighbor)
                print 'N pairs below threshold'
            else:  
                I1 = numpy.asarray(I1) # Model
                I2 = numpy.asarray(I2) # Lattice
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
                result = scipy.optimize.minimize(toBeMinimized_square, 
                                                 xStart, 
                                                 args = fixedParas, 
                                                 method = "Powell")
                print result.message
                
                ### RESULTS ###
                refinedScale  = result.x[0]
                refinedT_para = result.x[1]
                refinedT_perp = result.x[2]
                print refinedScale, refinedT_para, refinedT_perp
                
                scaled_I2 = (refinedScale * I2 * 
                             numpy.exp(-refinedT_para*qs_2D**2) * 
                             numpy.exp(-refinedT_perp*qs_rod**2))
                CC = correlate.Correlate(I1, scaled_I2)
                print 'N points: %d'%len(I1)
                print 'CC: %.4f'%CC
                
                delta_before = ((I1-I2)**2).sum()
                delta_after  = ((I1-scaled_I2)**2).sum()
                print 'DELTA FUNCTIONS: BEFORE %.2f AFTER %.2f'%(delta_before,
                                                                 delta_after)
                ###############################################################
                if CC >= CC_threshold:
                    nScaled = nScaled + 1

                    fOpen.write('L %s, K (L to M): %.2f, T (L to M): %.2f %.2f\n'
                                %(firstNeighbor, 
                                  refinedScale, 
                                  refinedT_para, 
                                  refinedT_perp))
                    print ('L %s, K (L to M): %.2f, T (L to M): %.2f %.2f\n'
                            %(firstNeighbor, 
                              refinedScale, 
                              refinedT_para, 
                              refinedT_perp))
                    
                    figureFlag = 0
                    if figureFlag == 1:
                        qs_rod_bins = numpy.arange(-1.5, 1.5, 0.01)
                        for i in range(0, len(qs_rod_bins)-1):
                            left = qs_rod_bins[i]
                            right = qs_rod_bins[i+1]
                            qs_2D_bin     = [qs_2D[j] 
                                             for j in range(0, len(qs_2D)) 
                                             if left <= qs_rod[j] < right]
                            I1_bin        = [I1[j]    
                                             for j in range(0, len(qs_2D)) 
                                             if left <= qs_rod[j] < right]
                            I2_bin        = [I2[j]    
                                             for j in range(0, len(qs_2D)) 
                                             if left <= qs_rod[j] < right]
                            scaled_I2_bin = [scaled_I2[j]    
                                             for j in range(0, len(qs_2D)) 
                                             if left <= qs_rod[j] < right]
                            if len(qs_2D_bin)>0:
                                matplotlib.pyplot.figure()
                                matplotlib.pyplot.title('qRod = %.2f - %.2f'
                                                        %(left, right))
                                matplotlib.pyplot.scatter(qs_2D_bin, 
                                                          I1_bin, 
                                                          marker= 'x', 
                                                          color='r', 
                                                          s=30)
                                matplotlib.pyplot.scatter(qs_2D_bin, 
                                                          I2_bin, 
                                                          marker= 'x', 
                                                          color='c', 
                                                          s=30)
                                matplotlib.pyplot.scatter(qs_2D_bin, 
                                                          scaled_I2_bin, 
                                                          marker= 'x', 
                                                          color='b', 
                                                          s=30)
                                matplotlib.pyplot.savefig('%s/lattice_%s_qRod_bin_%d.png'
                                                           %(figuresFolder,
                                                             firstNeighbor,
                                                             i),
                                                             dpi = 4*96)
                                matplotlib.pyplot.close()
                else:
                    refinedScale  = numpy.nan
                    refinedT_para = numpy.nan
                    refinedT_perp = numpy.nan
                    fOpen.write('L %s, CC below threshold.\n'%firstNeighbor) 
                    print 'L %s, CC below threshold.\n'%firstNeighbor
               
        LtoModel_vector.append([refinedScale, refinedT_para, refinedT_perp])   # L to M
    scaledFraction = float(nScaled)/nLattices
    fOpen.write('\nFraction of scaled lattices: %d/%d = %.2f.'%(nScaled,
                                                                nLattices,
                                                                scaledFraction)
                                                                              )  
    runTime = time.time() - startTime
    fOpen.write('\nElapsed time: %.1f s'%runTime)
    fOpen.close
    print 'Fraction of scaled lattices: %d/%d = %.2f.'%(nScaled,
                                                        nLattices,
                                                        scaledFraction)
                                                                              
   
    outputDir = '%s/r%s-scalesVsModel_anisotropic'%(inputFolder, runNumber)
    if not os.path.exists(outputDir):     
        os.mkdir(outputDir)           
    joblib.dump(LtoModel_vector, '%s/r%s-scalesVsModel_anisotropic.jbl'
                                  %(outputDir, runNumber))
           
    
if __name__ == "__main__":
    print "\n**** CALLING model_scaleVsModel_anisotropic ****"
    scalingFunction(sys.argv[1:])    