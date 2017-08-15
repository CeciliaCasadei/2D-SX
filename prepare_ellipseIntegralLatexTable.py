# -*- coding: utf-8 -*-
import sys
import getopt

def latexTable(myArguments):
    # READ INPUTS    
    try:
        optionPairs, leftOver = getopt.getopt(myArguments, "h", ["selectedRun="])
    except getopt.GetoptError:
        print 'Usage: python prepare_ellipseIntegralLatexTable.py --selectedRun <>selectedRun'
        sys.exit(2)   
    for option, value in optionPairs:
        if option == '-h':
            print 'Usage: python prepare_ellipseIntegralLatexTable.py --selectedRun <>selectedRun'
            sys.exit()
        elif option == "--selectedRun":
            runNumber = value.zfill(4)
        
    intensityValues = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/h_k_I_sum_circle_I_gauss_fixed_sigmas_I_sum_ellipse_x0_y0.txt'%runNumber, 'r')
    hs = []
    ks = []
    Is = []
    for line in intensityValues:
        h = int(line.split()[0])
        k = int(line.split()[1])
        I_ellipse = float(line.split()[4])
        hs.append(h)
        ks.append(k)
        Is.append(I_ellipse)
        print h, k, I_ellipse
    intensityValues.close()
    
    # SORTING
    # ORDER PEAKS ACCORDING TO INCREASING h
    for i in range(0, len(Is)):
        for j in range(i+1, len(Is)):
            if hs[i] > hs[j] or (hs[i] == hs[j] and ks[i] > ks[j]):
                h = hs[i]
                k = ks[i]
                I = Is[i]
                hs[i] = hs[j]
                ks[i] = ks[j]
                Is[i] = Is[j]
                hs[j] = h
                ks[j] = k
                Is[j] = I
                
    intensityValues_latexTable = open('./Output_r%s/Output_imageSums_moduleDisplacements_sigmaFits/h_k_I_sum_ellipse_latexTable.txt'%runNumber, 'w')
    
    for n in range(0, len(Is)):
        intensityValues_latexTable.write('%d & %d & %.2f \\\ \n'%(hs[n], ks[n], Is[n]  ))
    
    intensityValues_latexTable.close()

if __name__ == "__main__":
    print "\n**** CALLING prepare_ellipseIntegralLatexTable ****"
    latexTable(sys.argv[1:])   