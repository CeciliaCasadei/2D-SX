# -*- coding: utf-8 -*-
import getopt
import joblib
import random
import time
import numpy
import sys
import os

import scaling_calculateScaleFactor
import scaling_calculateScaleFactor_model

def scalingFunction(myArguments):
    
    # DEFAULTS
    runNumber = '' 
    deltaQrodThreshold = 0.001
    n_minThreshold = 8
    nTriangles = 100
    productThreshold = 0.25
    
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

    outputFolder = './Output_runMerging/transformAndScaleToModel_r%s'%runNumber
    
    # LOAD LATTICES LIST OF MATRICES: n h k qRod I Icorrected h_transformed k_transformed      
    myList = joblib.load('%s/spotsMatricesList-Transformed-r%s/r%s_transformedSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    nLattices = len(myList)
    
    # LOAD MODEL: h k qRod I        
    lattice_model = joblib.load('./Output_runMerging/startModel/lattice_model.jbl')
    lattice_model = numpy.asarray(lattice_model, dtype=numpy.float32)
   
    # SCALE LATTICES WITH RESPECT TO MODEL  
    startTime = time.time()
    #LtoModel_vector = []
    scaledLattices_list = []
    nScaled = 0
    
    fOpen = open('%s/r%s_scaling.txt'%(outputFolder, runNumber), 'w')
    fOpen.write('Total: %s lattices\n'%nLattices)
    fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
    fOpen.write('Min pairs n: %d\n'%n_minThreshold)
    fOpen.write('Product threshold: %f\n'%productThreshold)
    
    for firstNeighbor in range(0, nLattices):
        spots1stN = myList[firstNeighbor]  
        if spots1stN.shape[1] == 6:
            scale = numpy.nan
            fOpen.write('Lattice %s, Bad 1st N (non oriented)\n'%firstNeighbor)
            #LtoModel_vector.append(scale)
            continue
        
        print 'New 1st neighbor. Lattice %d'%firstNeighbor                
        nGood = 0
        nBad = 0
        
        spots1stN = numpy.asarray(spots1stN, dtype=numpy.float32) # n h k qRod I Icorrected h_transformed k_transformed
        n_min, scale_modelTo1stN = scaling_calculateScaleFactor_model.calculateScaleFactorFunction(lattice_model, spots1stN, deltaQrodThreshold)
        if n_min < n_minThreshold:
            scale = numpy.nan
            fOpen.write('Lattice %s, Bad 1st N\n'%firstNeighbor)
        else:
            secondShell = random.sample(range(nLattices), nTriangles)
            for secondNeighbor in secondShell:                    
                spots2ndN = myList[secondNeighbor]
                if spots2ndN.shape[1] == 8:
                    spots2ndN = numpy.asarray(spots2ndN, dtype=numpy.float32)
                    n_min, scale_1stNto2ndN = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots1stN, spots2ndN, 3*deltaQrodThreshold)
                    if n_min >= n_minThreshold:
                        n_min, scale_modelTo2ndN = scaling_calculateScaleFactor_model.calculateScaleFactorFunction(lattice_model, spots2ndN, deltaQrodThreshold)
                        scale_2ndNtoModel = float(1)/scale_modelTo2ndN
                        #n_min, scale_2ndNtoModel = scaling_calculateScaleFactor.calculateScaleFactorFunction(spots2ndN, lattice_modelForScaling, deltaQrodThreshold)
                        if n_min >= n_minThreshold:
                            product = scale_modelTo1stN*scale_1stNto2ndN*scale_2ndNtoModel
                            print product
                            if abs(product-1) <= productThreshold:
                                nGood = nGood + 1
                            else:
                                nBad = nBad + 1
                            
            if nGood+nBad >= 10 and float(nGood)/(nGood+nBad) >= 0.70:
                nScaled = nScaled + 1
                scale = float(1)/scale_modelTo1stN
                fOpen.write('Lattice %s, Scale: %.3f (nGood = %d, nBad = %d)\n'%(firstNeighbor, scale, nGood, nBad))
                spots1stN_scaled = []
                for spot in spots1stN:
                    scaledSpot = [spot[6], spot[7], spot[3], scale*spot[5]]   # h_transformed k_transformed qRod I_scaled
                    spots1stN_scaled.append(scaledSpot)
                spots1stN_scaled = numpy.asarray(spots1stN_scaled, dtype=numpy.float32)
            else:
                scale = numpy.nan
                fOpen.write('Lattice %s, Scaling: N\A (nGood = %d, nBad = %d)\n'%(firstNeighbor, nGood, nBad)) 
               
        #LtoModel_vector.append(scale) # lattice to model
        scaledLattices_list.append(spots1stN_scaled) # list of numpy 2D arrays: h_transformed k_transformed qRod I_scaled
    
    scaledFraction = float(nScaled)/nLattices
    fOpen.write('\nFraction of scaled lattices: %.2f.'%scaledFraction)    
    runTime = time.time() - startTime
    fOpen.write('\nIt took: %.1f s'%runTime)
    fOpen.close
                
    #joblib.dump(LtoModel_vector, '%s/r%s-scaling.jbl'%(outputFolder, runNumber))
    if not os.path.exists('%s/r%s_scaledLatticesList'%(outputFolder, runNumber)):
        os.mkdir('%s/r%s_scaledLatticesList'%(outputFolder, runNumber))
    joblib.dump(scaledLattices_list, '%s/r%s_scaledLatticesList/r%s_scaledLatticesList.jbl'%(outputFolder, runNumber, runNumber))
    
if __name__ == "__main__":
    print "\n**** CALLING model_scaleVsModel ****"
    scalingFunction(sys.argv[1:])    