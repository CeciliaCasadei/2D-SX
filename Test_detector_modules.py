# -*- coding: utf-8 -*-
import h5py
import numpy
import matplotlib.pyplot

import detectorModules

geometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5', 'r')
xGeometry = geometryFile['/x']   ### float32 ###
xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
yGeometry = geometryFile['/y']   ### float32 ###
yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)

modules = numpy.zeros((xGeometry_np.shape[0], xGeometry_np.shape[1]), dtype=numpy.int)
print xGeometry_np.shape
print yGeometry_np.shape
print modules.shape
modules = detectorModules.discontinuityMask(xGeometry_np, yGeometry_np, modules)
matplotlib.pyplot.imshow(modules, cmap='Greys')
matplotlib.pyplot.savefig('./modules.png', dpi=4*96)
matplotlib.pyplot.close()

row_interface = []
row_interface.append(0)
for i in range(0, modules.shape[0]):
    if modules[i, 0] == 1:
        row_interface.append(i)
row_interface.append(modules.shape[0]-1)

column_interface = []
column_interface.append(0)
for j in range(0, modules.shape[1]):
    if modules[0, j] == 1:
        column_interface.append(j)
column_interface.append(modules.shape[1]-1)

print row_interface
print column_interface