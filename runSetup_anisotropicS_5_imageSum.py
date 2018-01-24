# -*- coding: utf-8 -*-
import os
import sys


# SETUP FOR CURRENT RUN                                                                                                                                                                       
if len(sys.argv) != 2:
    print("[USAGE] %s tiltAngle" % sys.argv[0])
    sys.exit(-1)

tiltAngle = '%s'%(sys.argv[1])



# IMAGE SUM
binStep = 0.01             # A-1, in qz
binSize_factor = 3
resolutionLimit = 6.0      # A (2D)
halfWidth = 15             # pxls
imagesDirectoryName = '/afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS'
geometryFile = '/afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5'
nominalCell = 62.45        # A
wavelength = 1.48          # A
detectorDistance = 0.285   # m

flag = 0
if flag == 1:
    os.system('python imageSumming_tilted.py --binStep %f \
                                             --binSize_factor %f \
                                             --resolutionLimit %f \
                                             --halfWidth %d \
                                             --imagesDirectoryName %s \
                                             --geometryFile %s \
                                             --nominalCell %f \
                                             --wavelength %f \
                                             --tiltAngle %s \
                                             --detectorDistance %f'
                                             %(binStep, 
                                               binSize_factor, 
                                               resolutionLimit, 
                                               halfWidth, 
                                               imagesDirectoryName, 
                                               geometryFile, 
                                               nominalCell,
                                               wavelength,
                                               tiltAngle,
                                               detectorDistance))
    
flag = 0
if flag == 1:
    os.system('python imageSumming_tilted_spotWidth.py --halfWidth %d \
                                                       --tiltAngle %s'
                                                       %(halfWidth, 
                                                         tiltAngle))
                                                         
                                                         
flag = 0
if flag == 1:
    os.system('python imageSumming_tilted_spotWidthModeling.py --halfWidth %d \
                                                               --tiltAngle %s \
                                                               --nominalCell %f'
                                                               %(halfWidth, 
                                                                 tiltAngle,
                                                                 nominalCell))
                                                                 
 

ellipse_multiplicative_factor = 2.5
resolutionLimit = 6.0
thickness = 45
damping = 80
                                                                
flag = 1
if flag == 1:
    os.system('python imageSumming_tilted_finalIntegration.py --ellipse_multiplicative_factor %f \
                                                              --tiltAngle %s \
                                                              --resolutionLimit %f \
                                                              --nominalCell %f \
                                                              --thickness %f \
                                                              --damping %f'
                                                              %(ellipse_multiplicative_factor, 
                                                                tiltAngle,
                                                                resolutionLimit,
                                                                nominalCell,
                                                                thickness,
                                                                damping))
                                                                
