# -*- coding: utf-8 -*-
import h5py
import numpy

geometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5', 'r')
xGeometry = geometryFile['/x']   ### float32 ###
xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
yGeometry = geometryFile['/y']   ### float32 ###
yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
print 'Old'
print xGeometry_np.shape
print yGeometry_np.shape
print xGeometry_np[500:505, 700:705]

newGeometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r0201-images/geometry.h5', 'r')
new_xGeometry = newGeometryFile['/x']   ### float32 ###
new_xGeometry_np = numpy.asarray(new_xGeometry, dtype=numpy.float32)
new_yGeometry = newGeometryFile['/y']   ### float32 ###
new_yGeometry_np = numpy.asarray(new_yGeometry, dtype=numpy.float32)
print 'New'
print new_xGeometry_np.shape
print new_yGeometry_np.shape
print new_xGeometry_np[500:505, 700:705]

if numpy.array_equal(xGeometry_np, new_xGeometry_np):
    print 'New xGeo = Old xGeo'
else:
    print 'New xGeo differs from Old xGeo'
    
if numpy.array_equal(xGeometry_np, new_yGeometry_np):
    print 'New yGeo = Old xGeo'
else:
    print 'New yGeo differs from Old xGeo'
    
if numpy.array_equal(yGeometry_np, new_yGeometry_np):
    print 'New yGeo = Old yGeo'
else:
    print 'New yGeo differs from Old yGeo'