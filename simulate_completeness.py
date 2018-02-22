# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

import buildReciprocalLattice
import shannonSamplings

# PARAMETERS
cellSize = 62.45
wavelength = 1.48 #A
d = 45
resolutionLimit_3D = 5.3#phenix, all data 5.373
cutoff_2D = 6.0

# STANDARD PARAMETERS
tiltAngles_deg = [5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60.]
sinAlphas = numpy.linspace(-1, +1, 100, endpoint=True)
hmax = 20
kmax = 20

waveVector = 2 * numpy.pi / wavelength
c_star = (2*numpy.pi)/(2*d)
qLimit_3D = (2*numpy.pi)/resolutionLimit_3D
resolutionLimit_2D = resolutionLimit_3D

reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(cellSize,
                                                                          hmax,
                                                                          kmax,
                                                                          resolutionLimit_2D)

completenesses = []
matplotlib.pyplot.figure(figsize=(10, 15))
for tiltAngle_deg in tiltAngles_deg:
    tiltAngle = tiltAngle_deg * numpy.pi / 180 
    
    N_2D_exp_tot = 0
    N_sphere_tot = 0                                                                          
    for spot in reciprocalLattice:
        q_2D = float(spot[5])
        
        # 2D CRYSTAL    
        if q_2D > (2*numpy.pi)/cutoff_2D: 
            N_2D_exp = 0
        else:
            qs_rod = []
            for sinAlpha in sinAlphas:
                r = (numpy.cos(tiltAngle)**2 - 
                     (q_2D/waveVector)**2 - 
                      2*q_2D/waveVector*numpy.sin(tiltAngle)*sinAlpha)
                if r >= 0:
                    qrod = (
                             waveVector*
                             (
                              numpy.cos(tiltAngle)-
                              numpy.sqrt(r)
                              )
                            )
                    if abs(qrod) <= numpy.sqrt(qLimit_3D**2 - q_2D**2):
                        qs_rod.append(abs(qrod))
                else:
                    pass
                    #print 'Negative value in sqrt'
            qRod = max(qs_rod)
            samplings = shannonSamplings.get_shannonSamplings(c_star, qRod)
            samplings = samplings[1:-1]
            N_2D_exp = len(samplings)       
                    
        N_2D_exp_tot = N_2D_exp_tot + N_2D_exp
            
        # RECIPROCAL SPACE SPHERE
        qRod = numpy.sqrt(qLimit_3D**2 - q_2D**2)
        samplings = shannonSamplings.get_shannonSamplings(c_star, qRod)
        samplings = samplings[1:-1]
        N_sphere = len(samplings)
        N_sphere_tot = N_sphere_tot + N_sphere
        
    completeness = float(N_2D_exp_tot)/N_sphere_tot
    completenesses.append(completeness)
    print tiltAngle_deg, completeness
            
matplotlib.pyplot.figure()    
matplotlib.pyplot.scatter(tiltAngles_deg, completenesses)
axes = matplotlib.pyplot.gca()
axes.tick_params(axis='x', labelsize=12)
axes.tick_params(axis='y', labelsize=12)
axes.set_xlabel(r'$\eta$ (degrees)', fontsize = 14) 
axes.set_ylabel(r'Completeness',     fontsize = 14) 
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('./simulate_completeness.png', dpi=4*96)