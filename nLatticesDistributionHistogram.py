# -*- coding: utf-8 -*-
import matplotlib.pyplot
runNumber = '0201'
latticesFile = open('./Output_r%s/transformAndScale/spotsMatricesList-r%s/list-r%s.txt'%(runNumber, runNumber, runNumber), 'r')
lattices = list(latticesFile)

n_imagesDict = {}
n_imagesDict['0127'] = 968
n_imagesDict['0195'] = 308
n_imagesDict['0196'] = 440
n_imagesDict['0197'] = 484
n_imagesDict['0198'] = 440
n_imagesDict['0199'] = 528
n_imagesDict['0200'] = 440
n_imagesDict['0201'] = 528
n_images = int(n_imagesDict['%s'%runNumber])
print n_images
n_0 = 0
n_1 = 0
n_2 = 0
n_3 = 0
n_4 = 0
n_5 = 0
n_6 = 0
n_7 = 0
n_8 = 0

histogram = []
for image_n in range(0, n_images):
    n_lattices = 0
    for lattice in lattices:
        if int(lattice[80:84]) == image_n:
            print int(lattice[80:84])
            n_lattices = n_lattices + 1
    histogram.append(n_lattices)
    if n_lattices == 0:
        n_0 = n_0 + 1
    elif n_lattices == 1:
        n_1 = n_1 + 1
    elif n_lattices == 2:
        n_2 = n_2 + 1
    elif n_lattices == 3:
        n_3 = n_3 + 1
    elif n_lattices == 4:
        n_4 = n_4 + 1
    elif n_lattices == 5:
        n_5 = n_5 + 1
    elif n_lattices == 6:
        n_6 = n_6 + 1
    elif n_lattices == 7:
        n_7 = n_7 + 1
    elif n_lattices == 8:
        n_8 = n_8 + 1
        
print n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7, n_8
print 'n_images_tot', n_0+n_1+n_2+n_3+n_4+n_5+n_6+n_7+n_8
print 'n_lattices_tot', 0*n_0+1*n_1+2*n_2+3*n_3+4*n_4+5*n_5+6*n_6+7*n_7+8*n_8

matplotlib.pyplot.hist(histogram, bins=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
matplotlib.pyplot.gca().set_xlabel(r"Number of lattices", fontsize = 10, rotation = 'horizontal')
matplotlib.pyplot.gca().set_ylabel(r"Number of images",   fontsize = 10, rotation = 'vertical', labelpad = 3)
matplotlib.pyplot.savefig('./Output_r%s/N_lattices_after_processing_histogram.png'%runNumber)
matplotlib.pyplot.savefig('./Output_r%s/N_lattices_after_processing_histogram.pdf'%runNumber)

n0_perc = float(n_0)/n_images
n1_perc = float(n_1)/n_images
n2_perc = float(n_2)/n_images
n3_perc = float(n_3)/n_images
n4_perc = float(n_4)/n_images
n5_perc = float(n_5)/n_images
n6_perc = float(n_6)/n_images
print n0_perc, n1_perc, n2_perc, n3_perc, n4_perc, n5_perc, n6_perc
print n0_perc+n1_perc+n2_perc+n3_perc+n4_perc+n5_perc+n6_perc