# -*- coding: utf-8 -*-
# Count number of lattice-containing images, right after indexing.
import pickle

runNumber = '0127'
latticesDictionary_file = open('./Output_r%s/LatticeIndexing/r%s_allLatticesDictionary.pkl'%(runNumber, runNumber), 'rb')
latticesDictionary = pickle.load(latticesDictionary_file)
latticesDictionary_file.close()
nLattices = 0
imageNumbers = []
for i, j in latticesDictionary.items():
    nLattices = nLattices + 1
    imageNumbers.append(j.imageNumber)
print nLattices
print len(imageNumbers)
imageNumbers_unique = list(set(imageNumbers))
print len(imageNumbers_unique)