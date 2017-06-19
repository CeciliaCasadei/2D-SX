# -*- coding: utf-8 -*-
import h5py
import numpy
import matplotlib.pyplot

geometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5', 'r')
xGeometry = geometryFile['/x']   ### float32 ###
xGeometry_np = numpy.asarray(xGeometry, dtype=numpy.float32)
yGeometry = geometryFile['/y']   ### float32 ###
yGeometry_np = numpy.asarray(yGeometry, dtype=numpy.float32)
print 'Old'
print xGeometry_np.shape
print yGeometry_np.shape


newGeometryFile = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r0201-images/geometry.h5', 'r')
new_xGeometry = newGeometryFile['/x']   ### float32 ###
new_xGeometry_np = numpy.asarray(new_xGeometry, dtype=numpy.float32)
new_yGeometry = newGeometryFile['/y']   ### float32 ###
new_yGeometry_np = numpy.asarray(new_yGeometry, dtype=numpy.float32)
print 'New'
print new_xGeometry_np.shape
print new_yGeometry_np.shape

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
    

geometryFile_v1 = h5py.File('/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/2013-11-03_LCLS_CXI_LA52_geometry/cspad-cxia5213-v1.h5', 'r')
v1_xGeometry = geometryFile_v1['/x']   ### float32 ###
v1_xGeometry_np = numpy.asarray(v1_xGeometry, dtype=numpy.float32)
v1_yGeometry = geometryFile_v1['/y']   ### float32 ###
v1_yGeometry_np = numpy.asarray(v1_yGeometry, dtype=numpy.float32)
print 'Old'
print v1_xGeometry_np.shape
print v1_yGeometry_np.shape


if numpy.array_equal(xGeometry_np, v1_xGeometry_np):
    print 'v1_xGeo = Old xGeo'
else:
    print 'v1_xGeo differs from Old xGeo'
if numpy.array_equal(yGeometry_np, v1_yGeometry_np):
    print 'v1_yGeo = Old yGeo'
else:
    print 'v1_yGeo differs from Old yGeo'
    
matplotlib.pyplot.imshow(xGeometry_np, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./xGeo.png', dpi=4*96)
matplotlib.pyplot.close()
matplotlib.pyplot.imshow(yGeometry_np, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./yGeo.png', dpi=4*96)
matplotlib.pyplot.close()
matplotlib.pyplot.imshow(v1_xGeometry_np, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./xGeo_v1.png', dpi=4*96)
matplotlib.pyplot.close()
matplotlib.pyplot.imshow(v1_yGeometry_np, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./yGeo_v1.png', dpi=4*96)
matplotlib.pyplot.close()

xGeo_diff = xGeometry_np - v1_xGeometry_np
yGeo_diff = yGeometry_np - v1_yGeometry_np

matplotlib.pyplot.imshow(xGeo_diff, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./xGeo_diff.png')
matplotlib.pyplot.close()
matplotlib.pyplot.imshow(yGeo_diff, origin='lower')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.savefig('./yGeo_diff.png')
matplotlib.pyplot.close()

print '\n'
print xGeometry_np[500:505, 700:705]
print v1_xGeometry_np[500:505, 700:705] # diff 1.5 pxl
print xGeometry_np[500:505, 700:705] - v1_xGeometry_np[500:505, 700:705]
print yGeometry_np[500:505, 700:705] - v1_yGeometry_np[500:505, 700:705]
print '\n'
print xGeometry_np[1400:1405, 20:25]
print v1_xGeometry_np[1400:1405, 20:25] # diff 4 pxl
print xGeometry_np[1400:1405, 20:25] - v1_xGeometry_np[1400:1405, 20:25]
print yGeometry_np[1400:1405, 20:25] - v1_yGeometry_np[1400:1405, 20:25]
print '\n'
print xGeometry_np[1000:1005, 1000:1005]
print v1_xGeometry_np[1000:1005, 1000:1005] # diff 3 pxl
print xGeometry_np[1000:1005, 1000:1005] - v1_xGeometry_np[1000:1005, 1000:1005]