# -*- coding: utf-8 -*-
import sys
import getopt
import joblib

import makeOrbits
import simulate_resolution



def calculate_dataLimits_Function(myArguments):
    input_str_0 = '--inputFolder <inputFolder>'
    input_str_1 = '--resolutionLimit <resolutionLimit>'
    input_str_2 = '--cellSize <cellSize>'
    
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["inputFolder=", 
                                                                 "resolutionLimit=",
                                                                 "cellSize="])
    except getopt.GetoptError:
        print 'Usage: python dataLimits.py %s %s %s'%(input_str_0,
                                                      input_str_1,
                                                      input_str_2)
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python dataLimits.py %s %s %s'%(input_str_0,
                                                          input_str_1,
                                                          input_str_2)
            sys.exit()
        elif option == "--inputFolder":
            inputFolder = value
        elif option == "--resolutionLimit":
            resolutionLimit = float(value)
        elif option == "--cellSize":
            cellSize = float(value)
        
    
    # DEFINE ROD INDICES       
    orbits = makeOrbits.makeOrbitsFunction(resolutionLimit)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)
        
    print '%d Rods'%len(rodIndices)   
  
    ds = []
    N = 0
    for rod_hk in rodIndices:
        print rod_hk
        h = rod_hk[0]
        k = rod_hk[1]
        braggRodObject = joblib.load('%s/braggRodObjects/braggRodObject_%d_%d.jbl'
                                      %(inputFolder, h, k))
        
        experimental_q = braggRodObject.experimental_q
        for q in experimental_q:
            d = simulate_resolution.resolution(cellSize, h, k, q)
            ds.append(d)
     
    ds.sort()       
    dmin = min(ds)
    dmax = max(ds)
    N = len(ds)
    
    fOpen = open('%s/dataLimits.txt'%inputFolder, 'w')
    fOpen.write('d_max = %.4f A\n'%dmax)
    fOpen.write('d_min = %.4f A\n'%dmin)
    fOpen.write('N = %d \n'%N)
    
    nBins = 15
    binSize = int(N/nBins)
    for i in range(0, nBins+1):
        l = i*binSize
        r = (i+1)*binSize
        if r > len(ds):
            r = len(ds)
        ds_bin = ds[l : r]
        left = min(ds_bin)
        right = max(ds_bin)
        N_bin = len(ds_bin)
        fOpen.write('%10.4f %10.4f %d\n'%(left, right, N_bin))
        
    fOpen.close()

if __name__ == "__main__":
    print "\n**** CALLING calculate_dataLimits_Function ****"
    calculate_dataLimits_Function(sys.argv[1:])    
    