# -*- coding: utf-8 -*-
import getopt
import joblib
import time
import numpy
import sys
import os
import matplotlib
matplotlib.use('Agg') # Force matplotlib not to use any Xwindows backend.
import matplotlib.pyplot as m


import scaling_calculateScaleFactor
import correlate


def scalingFunction(myArguments):
    
    input_str_1 = '--runNumber <runNumber> --resolution_3D <resolution_3D>'
    input_str_2 = '--CC_threshold <CC_threshold> --model <model>'
    
    # DEFAULTS
    runNumber = '' 
    deltaQrodThreshold = 0.001
    n_minThreshold = 100
   
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, 
                                              "h", 
                                             ["runNumber=", 
                                              "resolution_3D=",
                                              "CC_threshold=",
                                              "model="])
    except getopt.GetoptError:
        print 'Usage: python model_scaleVsModel.py %s'%(input_str_1, 
                                                        input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_scaleVsModel.py %s'%(input_str_1,
                                                            input_str_2)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--resolution_3D":
            resolution_3D = float(value)
        elif option == "--CC_threshold":
            CC_threshold = float(value)
        elif option == "--model":
            model = value

    outputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%runNumber
    figuresFolder = '%s/scalingFigures'%outputFolder
    if not os.path.exists(figuresFolder):
        os.mkdir(figuresFolder)
    
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag i_unassembled j_unassembled
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'
                          %(outputFolder, runNumber, runNumber))
    nLattices = len(myList)
    print 'Total n of lattices: %d'%nLattices
    
    # LOAD MODEL: h k qRod I        
    lattice_model = joblib.load(model)
    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
   
    # SCALE LATTICES WITH RESPECT TO MODEL    
    LtoModel_vector = []
    nScaled = 0
    
    startTime = time.time()  
    fOpen = open('%s/r%s_scaling.txt'%(outputFolder, runNumber), 'w')
    fOpen.write('Total: %s lattices\n'%nLattices)
    fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
    fOpen.write('Min pairs n: %d\n'%n_minThreshold)
    fOpen.write('CC threshold: %f\n'%CC_threshold)
    
    for L in range(0, nLattices):
        spotsL = myList[L]  
        if spotsL[0, 4] == 0:
            scale = numpy.nan
            print 'L %s: Non oriented\n'%L
            fOpen.write('L %s: Non oriented\n'%L)            
        else:        
            print 'L %d'%L                

            n_pairs, \
            scale_modelToL, \
            I1,  \
            I2 = \
            scaling_calculateScaleFactor.calculateScaleFactorFunction(lattice_model, 
                                                                      spotsL, 
                                                                      deltaQrodThreshold,
                                                                      resolution_3D)
            if n_pairs < n_minThreshold:
                scale = numpy.nan
                fOpen.write('L %s: n_pairs model-lattice below threshold.\n'
                             %L)
            else:               
                I1 = numpy.asarray(I1) # MODEL
                I2 = numpy.asarray(I2) # LATTICE
                I1 = I1.flatten()
                I2 = I2.flatten()
                CC = correlate.Correlate(I1, I2)
                print 'N points: %d'%len(I1)
                print 'CC: %.4f'%CC
                
                if CC >= CC_threshold:
                    nScaled = nScaled + 1
                    scale = float(1)/scale_modelToL # LATTICE TO MODEL
                    fOpen.write('L %s, Scale (lattice to model): %.3f\n'
                                 %(L, scale))
                    print ('L %s, Scale (lattice to model): %.3f\n'
                            %(L, scale))
                    
                    figureFlag = 0
                    if figureFlag == 1:
                        m.figure()
                        m.scatter(I1, I2, color='b', s=2)
                        x = numpy.linspace(1.1*min(I1), 
                                           1.1*max(I1))
                        m.plot(x, scale_modelToL*x, 'r-')
                        m.axhline(y=0)
                        m.axvline(x=0)
    
                        m.title('CC = %.5f'%CC)
                        m.savefig('%s/scaleToModel_lattice%d.png'%(figuresFolder, 
                                                                   L))
                        m.close()    
                else:
                    scale = numpy.nan
                    fOpen.write('L %s, Scale: n/a\n'%L) 
                    print 'L %s, Scale: n/a\n'%L
               
        LtoModel_vector.append(scale) # lattice to model
        
    scaledFraction = float(nScaled)/nLattices
    fOpen.write('\nFraction of scaled lattices: %d/%d = %.2f.'%(nScaled,
                                                                nLattices,
                                                                scaledFraction))    
    runTime = time.time() - startTime
    fOpen.write('\nIt took: %.1f s'%runTime)
    fOpen.close
    print 'Fraction of scaled lattices: %d/%d = %.2f.'%(nScaled,
                                                        nLattices,
                                                        scaledFraction)
     
    if not os.path.exists('%s/r%s-scales'%(outputFolder, runNumber)):
        os.mkdir('%s/r%s-scales'%(outputFolder, runNumber))           
    joblib.dump(LtoModel_vector, '%s/r%s-scales/r%s-scales.jbl'
                                  %(outputFolder, runNumber, runNumber))
    
if __name__ == "__main__":
    print "\n**** CALLING model_scaleVsModel ****"
    scalingFunction(sys.argv[1:])    