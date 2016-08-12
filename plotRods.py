# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot

def plotRodsFunction(myArguments):
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber="])
    except getopt.GetoptError:
        print 'Usage: python plotRods.py --runNumber <runNumber>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRods.py --runNumber <runNumber>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)

    outputFolder = './Output_r%s/transformAndScale'%runNumber
    
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag 
    myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    
    # DEFINE ROD INDICES
    rodIndices = [[1, 0], [1, 1], [2, 0], [1, 2], [2, 1], [3, 0], [2, 2], [1, 3], [3, 1], [4, 0], [2, 3], [3, 2], [1, 4], [4, 1],
                  [5, 0], [3, 3], [2, 4], [4, 2], [1, 5], [5, 1], [6, 0], [3, 4], [4, 3], [2, 5], [5, 2], [1, 6], [6, 1],
                  [4, 4], [3, 5], [5, 3], [7, 0], [2, 6], [6, 2], [1, 7], [7, 1]]      
    
    # FOR EACH ROD, COLLECT P3-RELATED POINTS AND FRIEDEL MATES              
    for indices in rodIndices:
        hRod = indices[0]
        kRod = indices[1]
        Irod_vector = []
        Qrod_vector = []
        
        for latticeMatrix in myList:   
            latticeMatrix = numpy.asarray(latticeMatrix)
            if latticeMatrix[0, 4] == 0:                                         # Check lattice flag
                continue
            for spot in latticeMatrix:
                h_transformed = spot[0]
                k_transformed = spot[1]
                if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                    Irod_vector.append(spot[3])
                    Qrod_vector.append(spot[2])
                    
                if (h_transformed == -hRod and k_transformed == -kRod) or (h_transformed == hRod+kRod and k_transformed == -hRod) or (h_transformed == -kRod and k_transformed == hRod+kRod):
                    Irod_vector.append(spot[3])
                    Qrod_vector.append(-spot[2])

        # PLOT ROD
        matplotlib.pyplot.scatter(Qrod_vector, Irod_vector, marker='.', color='b')
        if not os.path.exists('%s/Rods'%outputFolder):
            os.mkdir('%s/Rods'%outputFolder)
        matplotlib.pyplot.savefig('%s/Rods/r%s_rod_%d_%d.png'%(outputFolder, runNumber, hRod, kRod))
        matplotlib.pyplot.close()


#        ROD GAPS
#        Run0198: Rod gaps are due to nan Intensities (i.e. integration circles touching detector gaps)        
#        matplotlib.pyplot.figure(facecolor = 'w')
#        (n, bins, patches) = matplotlib.pyplot.hist(Qrod_vector, bins=60)
#        matplotlib.pyplot.savefig('%s/Rods/r%s_rod_%d_%d_QrodDistribution.png'%(outputFolder, runNumber, hRod, kRod))
#        matplotlib.pyplot.close()
#        if runNumber == '0198':
#            if hRod == 1 and kRod == 3:
#                for i in range(0, len(Qrod_vector)):
#                    if 0.019 < Qrod_vector[i] < 0.026 or -0.026 < Qrod_vector[i] < -0.019:
#                        print 'Found qRod = %f I = %f'%(Qrod_vector[i], Irod_vector[i])
                     
if __name__ == "__main__":
    print "\n**** CALLING plotRods ****"
    plotRodsFunction(sys.argv[1:])    