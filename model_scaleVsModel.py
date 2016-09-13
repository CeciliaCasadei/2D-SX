# -*- coding: utf-8 -*-
import getopt
import joblib
import time
import numpy
import sys
import os
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot


import scaling_calculateScaleFactor

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

def scalingFunction(myArguments):

    # DEFAULTS
    runNumber = '' 
    deltaQrodThreshold = 0.001
    n_minThreshold = 100
    CC_threshold = 0.92
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python model_scaleVsModel.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_scaleVsModel.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_runMergingVsModel/transformAndScaleToModel_r%s'%runNumber
    figuresFolder = '%s/scalingFigures'%outputFolder
    if not os.path.exists(figuresFolder):
        os.mkdir(figuresFolder)
    
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag i_unassembled j_unassembled
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    nLattices = len(myList)
    print 'Total n of lattices: %d'%nLattices
    
    # LOAD MODEL: h k qRod I        
    lattice_model = joblib.load('./Output_runMerging/model/lattice_model.jbl')
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
    
    for firstNeighbor in range(0, nLattices):
        spots1stN = myList[firstNeighbor]  
        if spots1stN[0, 4] == 0:
            scale = numpy.nan
            print 'Lattice %s: Non oriented\n'%firstNeighbor
            fOpen.write('Lattice %s: Non oriented\n'%firstNeighbor)            
        else:        
            print 'Lattice %d'%firstNeighbor                

            n_pairs, scale_modelTo1stN, I1, I2 = scaling_calculateScaleFactor.calculateScaleFactorFunction(lattice_model, spots1stN, deltaQrodThreshold)
            if n_pairs < n_minThreshold:
                scale = numpy.nan
                fOpen.write('Lattice %s: n_pairs model-lattice below threshold.\n'%firstNeighbor)
            else:               
                ##########################################################################
                I1 = numpy.asarray(I1)
                I2 = numpy.asarray(I2)
                I1 = I1.flatten()
                I2 = I2.flatten()
                CC = Correlate(I1, I2)
                print 'N points: %d'%len(I1)
                print 'CC: %.4f'%CC
                
                
                matplotlib.pyplot.figure()
                matplotlib.pyplot.plot(I1, I2, 'bo')
                x = numpy.linspace(0, max(I1), num=float(max(I1))/0.001, endpoint = True)
                matplotlib.pyplot.plot(x, scale_modelTo1stN*x, 'r-')
                matplotlib.pyplot.title('CC = %.5f'%CC)
                matplotlib.pyplot.savefig('%s/scaleToModel_lattice%d.png'%(figuresFolder, firstNeighbor))
                matplotlib.pyplot.close()                    
                ##########################################################################
                if CC >= CC_threshold:
                    nScaled = nScaled + 1
                    scale = float(1)/scale_modelTo1stN
                    fOpen.write('Lattice %s, Scale (lattice to model): %.3f\n'%(firstNeighbor, scale))
                    print 'Lattice %s, Scale (lattice to model): %.3f\n'%(firstNeighbor, scale)
                else:
                    scale = numpy.nan
                    fOpen.write('Lattice %s, Scale: n/a\n'%firstNeighbor) 
                    print 'Lattice %s, Scale: n/a\n'%firstNeighbor
               
        LtoModel_vector.append(scale) # lattice to model
        
    scaledFraction = float(nScaled)/nLattices
    fOpen.write('\nFraction of scaled lattices: %.2f.'%scaledFraction)    
    runTime = time.time() - startTime
    fOpen.write('\nIt took: %.1f s'%runTime)
    fOpen.close
    print 'Total n of lattices: %d'%len(LtoModel_vector)  
    print 'Fraction of scaled lattices: %.2f.'%scaledFraction
     
    if not os.path.exists('%s/r%s-scales'%(outputFolder, runNumber)):
        os.mkdir('%s/r%s-scales'%(outputFolder, runNumber))           
    joblib.dump(LtoModel_vector, '%s/r%s-scales/r%s-scales.jbl'%(outputFolder, runNumber, runNumber))
    
if __name__ == "__main__":
    print "\n**** CALLING model_scaleVsModel ****"
    scalingFunction(sys.argv[1:])    