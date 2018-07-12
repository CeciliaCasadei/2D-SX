# -*- coding: utf-8 -*-
import os
import numpy

path = '/mnt/das-gpfs/home/casadei_c/work/casadei/Backup_tilted_Nov2016'

runs = ['0195', '0196', '0197', '0198', '0199', '0200', '0201']

a_s  = []
xb_s = []
yb_s = []
d_s  = []
 
for run in runs:
    path_r = '%s/Output_r%s/UnassembledImageProcessing'%(path, run)
    os.system('ls %s > temp_r%s.log'%(path_r, run) ) 
    fileList = open('temp_r%s.log'%run, 'r')
    for fileName in fileList:
        if 'Pro' in fileName:
            print fileName
            filePath = '%s/%s'%(path_r, fileName)
            filePath = filePath[:-1]
            file_o = open(filePath, 'r')
            file_o_list = list(file_o)
            n = len(file_o_list)
            for i in range(n):
                line = file_o_list[i]
                if 'Iteration: 2' in line:
                    a = file_o_list[i+2].split()[2]
                    xb = float(file_o_list[i+3].split()[2])
                    yb = float(file_o_list[i+4].split()[2])
                    d = numpy.sqrt(xb**2 + yb**2)
                    print a, xb, yb, d
                    a_s.append(float(a))
                    xb_s.append(xb)
                    yb_s.append(yb)
                    d_s.append(d)
                    
print len(a_s), len(xb_s), len(yb_s), len(d_s)

a_s  = numpy.asarray(a_s)
xb_s = numpy.asarray(xb_s)
yb_s = numpy.asarray(yb_s)
d_s  = numpy.asarray(d_s) 

a_avg = numpy.average(a_s)
a_std = numpy.std(a_s)

xb_avg = numpy.average(xb_s)
xb_std = numpy.std(xb_s)

yb_avg = numpy.average(yb_s)
yb_std = numpy.std(yb_s)

d_avg = numpy.average(d_s)
d_std = numpy.std(d_s)

print a_avg,  a_std
print xb_avg, xb_std
print yb_avg, yb_std
print d_avg, d_std
