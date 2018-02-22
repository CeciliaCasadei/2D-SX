# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

matplotlib.pyplot.figure(figsize=(10, 15))
wavelength = 1.48 #A
waveVector = 2 * numpy.pi / wavelength
q_2D = numpy.linspace(0, 2, 200)

sinAlphas = [-1, +1]
tilt_angles = [20, 45, 60] # Above 45 degrees, need to check all sin alphas
for tiltAngle_deg in tilt_angles:
    if tiltAngle_deg == 20:
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
matplotlib.pyplot.xticks(numpy.arange(0, 1.6, 0.4))


matplotlib.pyplot.gca().set_xlabel(r'$q_{\mathrm{2D}} (\AA^{-1})$', fontsize = 80)
matplotlib.pyplot.gca().set_ylabel(r'$q_{\mathrm{rod}} (\AA^{-1})$', fontsize = 80)
matplotlib.pyplot.gca().tick_params(labelsize=46, pad=5)
circle = matplotlib.pyplot.Circle((0, 0), 2*numpy.pi/5.3, alpha = 0.1)
matplotlib.pyplot.gca().add_artist(circle)
matplotlib.pyplot.gca().set_aspect('equal', adjustable='box')
matplotlib.pyplot.gca().set_xlim([0, 1.3])
matplotlib.pyplot.gca().set_ylim([-1.3, +1.3])
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('./qRod_vs_q2D_tilt_20_45_60.png', dpi=4*96)