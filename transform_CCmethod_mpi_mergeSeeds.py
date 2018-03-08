# -*- coding: utf-8 -*-
import joblib
import sys
import getopt
import os

def main(myArguments): 
       
    #READ COMMAND LINE ARGUMENTS
    str_1 = '--runNumber <runNumber>' 
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_mpi_merge.py %s'%(str_1)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_mpi_merge.py %s'%(str_1)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        
            
    transformationFolder = './Output_r%s/transformAndScale'%runNumber
    
    os.system('ls %s/r%s_orientations_*.jbl > %s/list.txt'%(transformationFolder, 
                                                            runNumber, 
                                                            transformationFolder))
    f = open('%s/list.txt'%transformationFolder, 'r')
    myList = list(f)
    f.close()
    
    nSeeds = len(myList)
    print 'N seeds: %d'%nSeeds
    
    transformations = []    
    for transformationsFile in myList:
        transformationsFile = transformationsFile.strip('\n')   
        seedResults = joblib.load(transformationsFile)
        print len(seedResults), ' lattices'
        transformations.append(seedResults)  
        
    
    joblib.dump(transformations, 
                '%s/r%s_orientations_allSeeds.jbl'
                 %(transformationFolder, 
                   runNumber))
                                      
    
    
if __name__ == "__main__":
    print "\n**** MERGING transform_CCmethod_mpi RESULTS: ALL SEEDS ****"
    main(sys.argv[1:])