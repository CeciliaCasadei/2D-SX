# -*- coding: utf-8 -*-
import numpy
import h5py
import pickle

import detectorModules

# EXTRACT GEOMETRY
geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5'     # same for all runs
geometryData = h5py.File(geometryFile, 'r')
xGeometry = geometryData['/x']   ### float32 ###
xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
yGeometry = geometryData['/y']   ### float32 ###
yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)

# EXTRACT DETECTOR MODULES BOUNDARIES
modules = numpy.zeros((xGeometry_np.shape[0], xGeometry_np.shape[1]), dtype=numpy.int)
modules = detectorModules.discontinuityMask(xGeometry_np, yGeometry_np, modules)

row_interface = []
row_interface.append(0)
for i in range(0, modules.shape[0]):
    if modules[i, 0] == 1:
        row_interface.append(i)
row_interface.append(modules.shape[0]-1)     #[0, 184, 369, 554, 739, 924, 1109, 1294, 1479]

column_interface = []
column_interface.append(0)
for j in range(0, modules.shape[1]):
    if modules[0, j] == 1:
        column_interface.append(j)
column_interface.append(modules.shape[1]-1)  #[0, 193, 387, 581, 775, 969, 1163, 1357, 1551]

print xGeometry_np.shape, yGeometry_np.shape
print row_interface
print column_interface

# LOOP ON MODULES
summary = []
for module_row in range(0, len(row_interface)-1):
    for module_column in range(0, len(column_interface)-1):
        bottom_bound = row_interface[module_row]
        top_bound = row_interface[module_row+1]
        left_bound = column_interface[module_column]
        right_bound = column_interface[module_column+1]
        print '\nMODULE: ', bottom_bound, top_bound, '---', left_bound, right_bound
        
        x0s_file = open('./Output_imageSums_moduleDisplacements_optimize/module_%d_%d_x0s.pkl'%(top_bound, right_bound), 'rb')
        x0s = pickle.load(x0s_file)
        x0s_file.close()
        
        y0s_file = open('./Output_imageSums_moduleDisplacements_optimize/module_%d_%d_y0s.pkl'%(top_bound, right_bound), 'rb')
        y0s = pickle.load(y0s_file)
        y0s_file.close()
        
        if len(x0s) > 0:
            avg_x0 = numpy.average(x0s)
            avg_y0 = numpy.average(y0s)
            std_x0 = numpy.std(x0s)
            std_y0 = numpy.std(y0s)
        else:
            avg_x0 = numpy.nan
            avg_y0 = numpy.nan
            std_x0 = numpy.nan
            std_y0 = numpy.nan
            
        print 'N: ', len(x0s)    
        print '%.2f %.2f   %.2f %.2f'%(avg_x0, std_x0, avg_y0, std_y0)
        
        module_results = [bottom_bound, top_bound, left_bound, right_bound, len(x0s), avg_x0, std_x0, avg_y0, std_y0]
        summary.append(module_results)
        
summary = numpy.asarray(summary)
print summary
print summary.shape