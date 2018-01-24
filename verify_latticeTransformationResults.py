# -*- coding: utf-8 -*-
runNumber = '0201'

n_I = 0
n_i = 0
n_p = 0
n_ip = 0
finalTransformations = open('./Output_r%s/transformAndScale/r%s-finalOrientations.txt'%(runNumber, runNumber), 'r')
for finalTransformation in finalTransformations:
    transformation = finalTransformation.split()[-1]
    if transformation == '0':
        n_I = n_I + 1
    elif transformation == '1':
        n_i = n_i + 1
    elif transformation == '2':
        n_p = n_p + 1
    elif transformation == '3':
        n_ip = n_ip + 1    
finalTransformations.close()

print n_I, n_i, n_p, n_ip
print n_I+n_i
print n_p+n_ip
print n_I+n_i+n_p+n_ip
print float(n_p+n_ip)/(n_I+n_i+n_p+n_ip)