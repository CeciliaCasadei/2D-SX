# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

import makeOrbits

matplotlib.pyplot.figure(figsize=(10, 15))
wavelength = 1.48 #A
waveVector = 2 * numpy.pi / wavelength
q_2D = numpy.linspace(0, 2, 200)

sinAlphas = [-1, +1]
tilt_angles = [20]
for tiltAngle_deg in tilt_angles:
    if tiltAngle_deg == 20:
        color = 'r'
    elif tiltAngle_deg == 45:
        color = 'g'
    else:
        color = 'b'
    for sinAlpha in sinAlphas:
        tiltAngle = tiltAngle_deg * numpy.pi / 180        
        qRod = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q_2D/waveVector)**2 - 2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha))
        matplotlib.pyplot.plot(q_2D, qRod, c=color)
        matplotlib.pyplot.plot(q_2D, -qRod, '--', c=color)
matplotlib.pyplot.legend()

matplotlib.pyplot.gca().set_xlabel(r'$q_{\mathrm{2D}} (\AA^{-1})$', fontsize = 28)
matplotlib.pyplot.gca().set_ylabel(r'$q_{\mathrm{rod}} (\AA^{-1})$', fontsize = 28)
matplotlib.pyplot.gca().tick_params(labelsize=22)
circle = matplotlib.pyplot.Circle((0, 0), 2*numpy.pi/5.4, alpha = 0.1)
matplotlib.pyplot.gca().add_artist(circle)
matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
matplotlib.pyplot.gca().set_xlim([0, 1.5])
matplotlib.pyplot.gca().set_ylim([-1.5, +1.5])
matplotlib.pyplot.tight_layout()

d = 45
c_star = (2*numpy.pi)/(2*d)
cellSize = 62.45 #62.45
#directCell = cellSize * numpy.matrix([[1, numpy.cos(2*numpy.pi/3)],[0, numpy.sin(2*numpy.pi/3)]]) # A
#reciprocalCellRows = 2* numpy.pi * directCell.I                                                   # A^(-1)
#
#
#
#orbits = makeOrbits.makeOrbitsFunction(6.0)
#print len(orbits)
#for orbit in orbits:
#    label = orbit.label
#    h = label[0]
#    k = label[1]
#    reciprocalVector = [h, k]*reciprocalCellRows
#    q_x = reciprocalVector[0,0]         # A^(-1)
#    q_y = reciprocalVector[0,1]         # A^(-1)
#    q_2D = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
tiltAngle = (20.0/180)*numpy.pi
sinAlpha = +1
resolution_3D = 5.4
q_max_3D = 2*numpy.pi/resolution_3D

n_obs = 0
n_tot = 0

hmax = 25
kmax = 25
lmax = 25
directCell = numpy.matrix([[cellSize, cellSize*numpy.cos(2*numpy.pi/3), 0  ],
                           [0,        cellSize*numpy.sin(2*numpy.pi/3), 0  ], 
                           [0,        0                               , 2*d]]) # A
reciprocalCellRows = 2* numpy.pi * directCell.I  
print reciprocalCellRows                                                 # A^(-1)
reciprocalLattice = []
i = 0
for h in range(-hmax, +hmax+1):
    print h
    for k in range(-kmax, +kmax+1):
        for l in range(0, +lmax+1):
            reciprocalVector = [h, k, l]*reciprocalCellRows
            q_x = reciprocalVector[0,0]         # A^(-1)
            q_y = reciprocalVector[0,1]         # A^(-1)
            q_z = reciprocalVector[0,2]
            q_2D = numpy.sqrt(q_x**2 + q_y**2)     # A^(-1)
            q_3D = numpy.sqrt(q_x**2 + q_y**2 + q_z**2)     # A^(-1)
            if q_2D == 0 and q_3D <= q_max_3D:
                n_tot = n_tot+1
            elif q_2D != 0 and q_3D <= q_max_3D :
                n_tot = n_tot+1
                q_z_max= (
                          waveVector*
                         (
                          numpy.cos(tiltAngle)-
                          numpy.sqrt(
                                     numpy.cos(tiltAngle)**2 - 
                                     (q_2D/waveVector)**2 - 
                                     2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha
                                     )
                          )
                          )
                if 0 <= q_z < q_z_max and (2*numpy.pi/q_2D)>6.0 :
                    n_obs = n_obs+1
                    matplotlib.pyplot.scatter(q_2D, q_z, s=1)
print n_obs, n_tot
print float(n_obs-66)/n_tot

matplotlib.pyplot.savefig('./qRod_vs_q2D_tilt_20_45_sinAlpha_pm1.png')

    
    
    