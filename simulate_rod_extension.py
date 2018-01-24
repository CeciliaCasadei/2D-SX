# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

matplotlib.pyplot.figure(figsize=(10, 15))
wavelength = 1.48 #A
waveVector = 2 * numpy.pi / wavelength
q_2D = numpy.linspace(0, 2, 200)

sinAlphas = [-1, +1]
tilt_angles = [15, 45]
for tiltAngle_deg in tilt_angles:
    if tiltAngle_deg == 15:
        color = 'r'
    elif tiltAngle_deg == 45:
        color = 'b'
    else:
        color = 'g'
    for sinAlpha in sinAlphas:
        tiltAngle = tiltAngle_deg * numpy.pi / 180        
        qRod = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q_2D/waveVector)**2 - 2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha))
        if sinAlpha == +1:
            matplotlib.pyplot.plot(q_2D, qRod, c=color)
            matplotlib.pyplot.plot(q_2D, -qRod, c=color)
#        else:
#            matplotlib.pyplot.plot(q_2D, qRod, '--', c=color)
#            matplotlib.pyplot.plot(q_2D, -qRod, '--', c=color)
            
matplotlib.pyplot.legend()

matplotlib.pyplot.gca().set_xlabel(r'$q_{\mathrm{2D}} (\AA^{-1})$', fontsize = 40)
matplotlib.pyplot.gca().set_ylabel(r'$q_{\mathrm{rod}} (\AA^{-1})$', fontsize = 40)
matplotlib.pyplot.gca().tick_params(labelsize=32, pad=20)
circle = matplotlib.pyplot.Circle((0, 0), 2*numpy.pi/3.5, alpha = 0.1)
matplotlib.pyplot.gca().add_artist(circle)
matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
matplotlib.pyplot.gca().set_xlim([0, 2.0])
matplotlib.pyplot.gca().set_ylim([-2.0, +2.0])
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('./qRod_vs_q2D_tilt_15_45_sinAlpha_p1.png')

#tiltAngle = 35 * numpy.pi / 180     
#q_2D = 1.14
#sinAlpha = 1
#qRod_plus = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q_2D/waveVector)**2 - 2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha))
#q_3D = numpy.sqrt(q_2D**2 + qRod_plus**2)
#sinAlpha = -1
#qRod_minus = waveVector*(numpy.cos(tiltAngle)-numpy.sqrt(numpy.cos(tiltAngle)**2 - (q_2D/waveVector)**2 - 2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha))
#print q_2D, qRod_plus, qRod_minus, q_3D, 2*numpy.pi/q_3D