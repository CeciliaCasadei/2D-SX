# -*- coding: utf-8 -*-
runNumber = '0127'
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

intensityValues_latexTable.write('\\begin{center}\n')
intensityValues_latexTable.write('\\begin{tabular}{ |l|l|l| } \n')
intensityValues_latexTable.write('h & k & I \\\ \n\\hline\n')
#\begin{tabular}{ |c|c|c| } 
 #\hline
 #cell1 & cell2 & cell3 \\ 
 #cell4 & cell5 & cell6 \\ 
 
#\end{tabular}
#\end{center}




for n in range(0, len(Is)):
    intensityValues_latexTable.write('%d & %d & %.2f \\\ \n'%(hs[n], ks[n], Is[n]  ))


intensityValues_latexTable.write('\\end{tabular}\n')

intensityValues_latexTable.write('\\end{center}\n')
intensityValues_latexTable.close()