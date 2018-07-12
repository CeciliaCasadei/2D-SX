# -*- coding: utf-8 -*-
import joblib
import numpy 
import time
import sys
import getopt
import os

### CYTHON MODULES ###
import transform_calculateCCs

    
def main(myArguments): 
    identity = numpy.matrix([[1, 0],[0, 1]])
    
    # DEFAULTS:
    deltaQrodThreshold = 0.001
    n_minThreshold = 6
    
    str_input_1 = '--runNumber <runNumber> --model <model> --dQrod <dQrod>'
    str_input_2 = '--nMin <nMin> --nLattices <nLattices>'
    
    # READ COMMAND LINE ARGUMENTS
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=",
                                                               "model=", 
                                                               "dQrod=", 
                                                               "nMin=", 
                                                               "nLattices="])
    except getopt.GetoptError:
        print 'Usage: python model_transformVsModel.py %s %s %s'%(str_input_1, 
                                                                  str_input_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python model_transformVsModel.py %s %s %s'%(str_input_1, 
                                                                      str_input_2)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--model":
            lattice_model = value   
        elif option == "--dQrod":
            deltaQrodThreshold = float(value)
        elif option == "--nMin":
            n_minThreshold = int(value)
        elif option == "--nLattices":
            nUsedLattices = int(value)
    
    newFolder = './Output_runMergingVsModel'
    transformationFolder = '%s/transformAndScaleToModel_r%s'%(newFolder, 
                                                              runNumber)                                                           
    if not os.path.exists(transformationFolder):
        os.mkdir(transformationFolder)

    # LOAD MODEL: h k qRod I
    lattice_model = joblib.load('%s'%lattice_model)  
            
    # LOAD LATTICES LIST OF MATRICES: h k qRod I flag=1 i j
    listPath = './Output_r%s/transformAndScale/spotsMatricesList-r%s'%(runNumber, 
                                                                       runNumber)
    myList = joblib.load('%s/r%s_spotsMatricesList.jbl'%(listPath, runNumber))
    nLattices = len(myList)
    
    if not 'nUsedLattices' in locals():        
        nUsedLattices = int(nLattices)
    print 'Number of lattices to transform: %d'%nUsedLattices
    
    startTime = time.time()
    
    # ORIENT LATTICES 0 - nUsedLattices-1 WITH RESPECT TO MODEL
    nOriented = 0
    LtoModel_vector = []
    
    fOpen = open('%s/r%s-orientationsVsModel.txt'%(transformationFolder, 
                                                   runNumber), 
                                                   'w')
    fOpen.write('Total: %s lattices\n'%nLattices)
    fOpen.write('Delta qRod: %f\n'%deltaQrodThreshold)
    fOpen.write('Min pairs n: %d\n'%n_minThreshold)
   
    for L in range(0, nUsedLattices):
        print '\nLattice %d'%L                
        
        spotsL = myList[L]   
        print spotsL.shape            
        n_min, \
        avg_CCs \
        = transform_calculateCCs.determineTransformation(lattice_model, 
                                                         spotsL, 
                                                         deltaQrodThreshold)
        if n_min < n_minThreshold:
            transformation_ModelToL = numpy.nan
            print 'n pairs below threshold for lattice ', L
        else:
            print 'n pairs above threshold (model to L): %d'%n_min
            transformation_ModelToL = avg_CCs.index(max(avg_CCs))
            nOriented = nOriented + 1
    
        fOpen.write('Lattice %s, Orientation: %s\n'%(L, 
                                                     transformation_ModelToL))
        LtoModel_vector.append(transformation_ModelToL)
        
    fractionOriented = float(nOriented)/nLattices
    fOpen.write('\nFraction of transformed lattices: %d/%d = %.2f.'
                 %(nOriented, nLattices, fractionOriented))
    runTime = time.time() - startTime
    fOpen.write('\nElapsed time: %.1f s'%runTime)
    fOpen.close
            
    joblib.dump(LtoModel_vector, 
                '%s/r%s_orientationsVsModel.jbl'%(transformationFolder, 
                                                  runNumber))
    
if __name__ == "__main__":
    print "\n**** CALLING model_transformVsModel ****"
    main(sys.argv[1:])