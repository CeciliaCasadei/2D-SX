# -*- coding: utf-8 -*-
import joblib
import sys
import getopt
import os

def main(myArguments): 
       
    #READ COMMAND LINE ARGUMENTS
    str_1 = '--runNumber <runNumber> --seed_id <seed_id>' 
    try:
        optionPairs, leftOver = getopt.getopt(myArguments,"h",["runNumber=", 
                                                               "seed_id="])
    except getopt.GetoptError:
        print 'Usage: python transform_CCmethod_mpi_merge.py %s'%(str_1)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python transform_CCmethod_mpi_merge.py %s'%(str_1)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--seed_id":
            seed_id = float(value)
        
            
    transformationFolder = './Output_r%s/transformAndScale'%runNumber
    
    os.system('ls %s/r%s_orientations_L_*.jbl > %s/list.txt'%(transformationFolder, 
                                                              runNumber, 
                                                              transformationFolder))
    f = open('%s/list.txt'%transformationFolder, 'r')
    myList = list(f)
    f.close()
    nLattices = len(myList)
    print 'N processed lattices: %d'%nLattices
    
    fOpen = open('%s/r%s_orientations_%d.txt'%(transformationFolder, 
                                               runNumber, 
                                               seed_id), 
                                               'w')
    
    transformations = []    
    for myLatticeFile in myList:
        myLatticeFile = myLatticeFile.strip('\n')   
        latticeResults = joblib.load(myLatticeFile)
        L = int(latticeResults[0])
        T = latticeResults[1]
        fOpen.write('%10d %s\n'%(L, T))
        transformations.append(T)  
        
    fOpen.close()
    
    joblib.dump(transformations, 
                '%s/r%s_orientations_%d.jbl'
                 %(transformationFolder, 
                   runNumber, 
                   seed_id))
                                      
    os.system('rm %s/r%s_orientations_L_*.jbl'%(transformationFolder, 
                                                runNumber))
    
    
    
if __name__ == "__main__":
    print "\n**** MERGING transform_CCmethod_mpi RESULTS ****"
    main(sys.argv[1:])