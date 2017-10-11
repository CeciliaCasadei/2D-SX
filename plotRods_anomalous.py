# -*- coding: utf-8 -*-
import sys
import getopt
import os
import joblib
import numpy
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot


import makeOrbits

def plotRodsFunction(myArguments):
    
    # DEFAULTS
    outputFolder = ''
    resolutionLimit = 7.0
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["runNumber=", "outputFolder="])
    except getopt.GetoptError:
        print 'Usage: python plotRods.py --runNumber <runNumber> --outputFolder <outputFolder>'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python plotRods.py --runNumber <runNumber> --outputFolder <outputFolder>'
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--outputFolder":
            outputFolder = value

    print runNumber
    if outputFolder == '':
        outputFolder = './Output_r%s/transformAndScale'%runNumber
    
        
    print 'OUTPUT FOLDER ', outputFolder
    outputFile = open('%s/mergedIntensities_h_k_qRod_l_avgI_sigmaI.txt'%outputFolder, 'w')
    
    # FOLDERS   
    if not os.path.exists('%s/Histograms'%outputFolder):
        os.mkdir('%s/Histograms'%outputFolder)
    
    # LOAD LATTICES LIST OF MATRICES: h_transformed k_transformed qRod Iscaled flag i_unassembled j_unassembled scale
    myList = joblib.load('%s/spotsMatricesList-Scaled-r%s/r%s_scaledSpotsMatricesList.jbl'%(outputFolder, runNumber, runNumber))
    
    # DEFINE ROD INDICES 
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    orbitIndices = []
    for orbit in orbits:
        orbit_label = orbit.label        
        orbitIndices.append(orbit_label)
        
    print '%d orbits'%len(orbitIndices) # (h, k), (-h, -k), (k, h), (-k, -h)
    

    #FOR EACH ORBIT, COLLECT P3-RELATED POINTS AND SEPARATE FRIEDEL MATES              
    for indices in orbitIndices:
        hRod = indices[0]
        kRod = indices[1]
        Irod_vector_posQ = []
        Qrod_vector_posQ = []
        Irod_vector_negQ = []
        Qrod_vector_negQ = []
        
        
        for latticeMatrix in myList:   
            latticeMatrix = numpy.asarray(latticeMatrix)
            if latticeMatrix[0, 4] == 0:                                         # Check lattice flag
                continue
            for spot in latticeMatrix:
                h_transformed = spot[0]
                k_transformed = spot[1]
                q_transformed = spot[2]
                if (h_transformed == hRod and k_transformed == kRod) or (h_transformed == -hRod-kRod and k_transformed == hRod) or (h_transformed == kRod and k_transformed == -hRod-kRod):
                    if q_transformed > 0:
                        Irod_vector_posQ.append(spot[3])
                        Qrod_vector_posQ.append(spot[2])
                    elif q_transformed < 0:
                        Irod_vector_negQ.append(spot[3])
                        Qrod_vector_negQ.append(spot[2])
                    else:
                        print 'ERROR'
                    
                
        
        
        # LOG RESULTS
        Irod_vector_posQ_cleaned = [Irod_vector_posQ[i] for i in range(0, len(Irod_vector_posQ)) if not numpy.isnan(Irod_vector_posQ[i])]
        Irod_vector_negQ_cleaned = [Irod_vector_negQ[i] for i in range(0, len(Irod_vector_negQ)) if not numpy.isnan(Irod_vector_negQ[i])]
        Qrod_vector_posQ_cleaned = [Qrod_vector_posQ[i] for i in range(0, len(Qrod_vector_posQ)) if not numpy.isnan(Irod_vector_posQ[i])]
        Qrod_vector_negQ_cleaned = [Qrod_vector_negQ[i] for i in range(0, len(Qrod_vector_negQ)) if not numpy.isnan(Irod_vector_negQ[i])]
        
        avg_Intensity_posQ = numpy.average(Irod_vector_posQ_cleaned)
        std_Intensity_posQ = numpy.std(Irod_vector_posQ_cleaned)
        avg_Intensity_negQ = numpy.average(Irod_vector_negQ_cleaned)
        std_Intensity_negQ = numpy.std(Irod_vector_negQ_cleaned)
        
        avg_Qrod_posQ = numpy.average(Qrod_vector_posQ_cleaned)
        std_Qrod_posQ = numpy.std(Qrod_vector_posQ_cleaned)
        avg_Qrod_negQ = numpy.average(Qrod_vector_negQ_cleaned)
        std_Qrod_negQ = numpy.std(Qrod_vector_negQ_cleaned)
        
        
        
        
        
        outputFile.write('\n%3d %3d %7.4f %7.2f %7.2f %7.2f %7.1f'%( hRod,  kRod, avg_Qrod_posQ, std_Qrod_posQ, avg_Intensity_posQ, std_Intensity_posQ, (avg_Intensity_posQ/std_Intensity_posQ)))
        outputFile.write('\n%3d %3d %7.4f %7.2f %7.2f %7.2f %7.1f'%( hRod,  kRod, avg_Qrod_negQ, std_Qrod_negQ, avg_Intensity_negQ, std_Intensity_negQ, (avg_Intensity_negQ/std_Intensity_negQ)))
         
        try:
            # PLOT ORBIT HISTOGRAMS
            Irod_vector_posQ_cleaned = numpy.ravel(Irod_vector_posQ_cleaned)
            Irod_vector_negQ_cleaned = numpy.ravel(Irod_vector_negQ_cleaned)
            
            matplotlib.pyplot.figure(facecolor = 'w')
            matplotlib.pyplot.hist(Irod_vector_posQ_cleaned)
            matplotlib.pyplot.xlabel(r"I$_{\rm scaled}$ (photons)")
            matplotlib.pyplot.ylabel(r"N$_{\rm measured}$")        
            matplotlib.pyplot.savefig('%s/Histograms/Histogram_r%s_spot_%d_%d_%.2f.png'%(outputFolder, runNumber, hRod, kRod, avg_Qrod_posQ), dpi = 300, facecolor = 'w')
            matplotlib.pyplot.close()
                   
            matplotlib.pyplot.figure(facecolor = 'w')
            matplotlib.pyplot.hist(Irod_vector_negQ_cleaned)
            matplotlib.pyplot.xlabel(r"I$_{\rm scaled}$ (photons)")
            matplotlib.pyplot.ylabel(r"N$_{\rm measured}$")
            matplotlib.pyplot.savefig('%s/Histograms/Histogram_r%s_spot_%d_%d_%.2f.png'%(outputFolder, runNumber, hRod, kRod, avg_Qrod_negQ), dpi = 300, facecolor = 'w')
            matplotlib.pyplot.close()
        except:
            print 'HISTOGRAM ERROR'
        
            
    outputFile.close()        

if __name__ == "__main__":
    print "\n**** CALLING plotRods ****"
    plotRodsFunction(sys.argv[1:])    