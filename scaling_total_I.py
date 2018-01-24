# -*- coding: utf-8 -*-
import getopt
import joblib
import numpy
import sys
import matplotlib.pyplot
import os

import simulate_resolution

def calculate_latticeAvgI(lattices, 
                          cellSize, 
                          resolution_3D):
    avg_Is = []
    i = -1
    for lattice in lattices:
        i = i+1
        I = numpy.nan
        # CHECK FLAG
        if lattice[0, 4] == 1:
            n = 0
            I = 0
            for spot in lattice:
                h = spot[0]
                k = spot[1]
                qRod = spot[2]
                resolution = simulate_resolution.resolution(cellSize, 
                                                            h, 
                                                            k, 
                                                            qRod)
                if resolution >= resolution_3D:
                    if not numpy.isnan(spot[3]):
                        n = n + 1
                        I = I + spot[3]
            I = I/n
            print 'LATTICE %d AVG I: %.2f ph, FROM %d SPOTS'%(i, I, n)
        avg_Is.append(I)
    return avg_Is
    
    
def calculate_datasetTotI(lattices, 
                          cellSize, 
                          resolution_3D):
    I = 0
    n = 0
    for lattice in lattices:
        # CHECK FLAG
        if lattice[0, 4] == 1:
            for spot in lattice:
                h = spot[0]
                k = spot[1]
                qRod = spot[2]
                resolution = simulate_resolution.resolution(cellSize, 
                                                            h, 
                                                            k, 
                                                            qRod)
                if resolution >= resolution_3D:
                    if not numpy.isnan(spot[3]):
                        I = I + spot[3]  
                        n = n+1
    return n, I

def plotDistrivbution(data, avg, title, xlabel, ylabel, figureName):
    matplotlib.pyplot.figure()
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.hist(data)
    matplotlib.pyplot.axvline(x=avg)
    matplotlib.pyplot.gca().set_xlabel(r'%s'%xlabel, 
                                       fontsize = 20, 
                                       rotation = 'horizontal')
    matplotlib.pyplot.gca().set_ylabel(r'%s'%ylabel, 
                                       fontsize = 20, 
                                       rotation = 'vertical')
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.close()

def scalingFunction(myArguments):    
        
    input_str_1 = '--runNumber <runNumber>'
    input_str_2 = '--resolution_3D <resolution_3D>'
    input_str_3 = '--cellSize <cellSize>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", 
                                             ["runNumber=", 
                                              "resolution_3D=",
                                              "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python scaling_total_I.py %s %s %s'%(input_str_1, 
                                                           input_str_2,
                                                           input_str_3)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python scaling_total_I.py %s %s %s'%(input_str_1, 
                                                               input_str_2,
                                                               input_str_3)
            sys.exit()
        elif option == "--runNumber":
            runNumber = value.zfill(4)
        elif option == "--resolution_3D":
            resolution_3D = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
                     
    outputFolder = './Output_r%s/transformAndScale'%runNumber        
    print 'OUTPUT FOLDER: ', outputFolder
       
    # LOAD LATTICES LIST OF MATRICES: h_t k_t qRod I flag i_unas j_unas
    folder = '%s/spotsMatricesList-Transformed-r%s'%(outputFolder, runNumber)
    myList = joblib.load('%s/r%s_transformedSpotsMatricesList.jbl'
                          %(folder, runNumber))
    print 'TOTAL N LATTICES: ', len(myList)
    
    # CALCULATE DATASET TOT LOW RES I
    nUsed, tot_I = calculate_datasetTotI(myList, cellSize, resolution_3D)
    print 'DATASET TOT I BEFORE SCALING: %.2f (%d spots)'%(tot_I, nUsed)
    
    # CALCULATE AVG LOW RES I PER LATTICE
    avg_Is = calculate_latticeAvgI(myList, cellSize, resolution_3D)
    
    # CALCULATE AVERAGE ON LATTICES
    avg_Is_cleaned = [avg_Is[j] for j in range(len(avg_Is)) 
                                if not numpy.isnan(avg_Is[j])]   
    I_avg_lattice = numpy.average(avg_Is_cleaned)  
    print 'N LATTICES TO SCALE: %d/%d '%(len(avg_Is_cleaned), len(avg_Is))
    print 'I AVG ON LATTICES: %.2f ph '%(I_avg_lattice)
    
    # PLOT AVG LOW RES I DISTRIBUTION  
    title = 'Distribution of lattice average low resolution I before scaling.'
    xlabel = "Average I (photons)"
    ylabel = "N of lattices"
    figureName = '%s/r%s_avg_lowRes_I.png'%(outputFolder, runNumber)
    plotDistrivbution(avg_Is_cleaned, 
                      I_avg_lattice, 
                      title, 
                      xlabel, 
                      ylabel, 
                      figureName)
    
    # DETERMINE SCALE FACTORS
    Ks = []
    for avg_I in avg_Is:
        if numpy.isnan(avg_I):
            K = numpy.nan
        else:
            K = I_avg_lattice/avg_I
        Ks.append(K)
    print 'N LATTICES: ', len(Ks)
    Ks_cleaned = [Ks[j] for j in range(len(Ks)) if not numpy.isnan(Ks[j])]  
    avg_K = numpy.average(Ks_cleaned)
    print 'Average K: %.2f'%avg_K
    
    # PLOT K DISTRIBUTION
    title = 'Distribution of scale constant K.'
    xlabel = "K"
    ylabel = "N of lattices"
    figureName = '%s/r%s_Ks.png'%(outputFolder, runNumber)
    plotDistrivbution(Ks_cleaned, 
                      avg_K,
                      title, 
                      xlabel, 
                      ylabel, 
                      figureName)
    
    # APPLY SCALE FACTORS
    scaledLattices = []
    nScaled = 0
    for i in range(len(myList)):
        lattice = myList[i]
        scale = Ks[i]
        if numpy.isnan(scale):
            unscaledLattice = []
            for spot in lattice:
                unscaledSpot = [spot[0],        # h
                                spot[1],        # k
                                spot[2],        # qRod
                                spot[3],        # I, unscaled
                                spot[4],        # flag = 0
                                spot[5],        # i_unassembled
                                spot[6],        # j_unassembled
                                scale]          # nan
                unscaledLattice.append(unscaledSpot)
            unscaledLattice = numpy.asarray(unscaledLattice)
            scaledLattices.append(unscaledLattice)
        else:
            nScaled = nScaled + 1
            scaledLattice = []
            for spot in lattice:
                scaledSpot = [spot[0],        # h
                              spot[1],        # k
                              spot[2],        # qRod
                              scale*spot[3],  # I
                              spot[4],        # flag
                              spot[5],        # i_unassembled
                              spot[6],        # j_unassembled
                              scale]          # K
                scaledLattice.append(scaledSpot)
            scaledLattice = numpy.asarray(scaledLattice)
            scaledLattices.append(scaledLattice)
    print 'N LATTICES: ', len(scaledLattices)
    print 'N SCALED: ', nScaled
    
    # SAVE SCALED LATTICES
    outputPath = '%s/spotsMatricesList-Scaled-r%s'%(outputFolder, runNumber)
    if not os.path.exists(outputPath):
        os.mkdir(outputPath)
    joblib.dump(scaledLattices, 
                '%s/r%s_scaledSpotsMatricesList.jbl'%(outputPath, runNumber))
   
    # CALCULATE AVG LOW RES I PER LATTICE, AFTER SCALING      
    avg_Is_afterS = calculate_latticeAvgI(scaledLattices, 
                                          cellSize, 
                                          resolution_3D)    
    avg_Is_afterS_cleaned = [avg_Is_afterS[j] 
                             for j in range(len(avg_Is_afterS)) 
                             if not numpy.isnan(avg_Is_afterS[j])]   
    I_avg_lattice_afterS = numpy.average(avg_Is_afterS_cleaned)  
    print 'I AVG ON LATTICES: %.2f ph (AFTER SCALING)'%(I_avg_lattice_afterS)
    
    # CALCULATE DATASET TOT LOW RES I
    nUsed_afterS, tot_I_afterS = calculate_datasetTotI(scaledLattices, 
                                                       cellSize, 
                                                       resolution_3D)
    print 'DATASET TOT I AFTER SCALING: %.2f (%d spots)'%(tot_I_afterS, 
                                                          nUsed_afterS)

    

if __name__ == "__main__":
    print "\n**** CALLING scaling ****"
    scalingFunction(sys.argv[1:])    