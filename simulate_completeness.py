# -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot

import buildReciprocalLattice
import shannonSamplings

from bins import binLimits

def get_qRod(sinAlphas, tiltAngle, q_2D, waveVector, qRodLimit):
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
            if abs(qrod) <= qRodLimit:
                qs_rod.append(abs(qrod))
    qRod = max(qs_rod)
    return qRod
    
def getN(samplingSpace, qRodMax, q_2D, bin_q3D_low, bin_q3D_high):
    samplings = shannonSamplings.get_shannonSamplings(samplingSpace, 
                                                      qRodMax)
    samplings = samplings[1:-1]
    
    N = 0
    for qRod_s in samplings:
        q_3D = numpy.sqrt(q_2D**2 + qRod_s**2)
        if bin_q3D_low <= q_3D < bin_q3D_high:
            N = N+1
        
    return N, len(samplings)
    
def getBinCompleteness(reciprocalLattice, 
                       qLimit_3D, 
                       sinAlphas, 
                       tiltAngle,
                       waveVector,
                       samplingSpace,
                       q3D_low, 
                       q3D_high):
    N_2D_exp_tot = 0
    N_sphere_tot = 0     
                                                                     
    for spot in reciprocalLattice:
        q_2D = float(spot[5])
        qRod_limit_resSphere = numpy.sqrt(qLimit_3D**2 - q_2D**2)
        
        # 2D CRYSTAL                  
        qRod = get_qRod(sinAlphas, 
                        tiltAngle, 
                        q_2D, 
                        waveVector, 
                        qRod_limit_resSphere)
#        print 'q2D = %6.2f qRod_sphere = %6.2f qRod_2Dxtal = %.2f'%(q_2D, 
#                                                                    qRod_limit_resSphere, 
#                                                                    qRod)
        N, N_samplings = getN(samplingSpace, 
                              qRod, 
                              q_2D, 
                              q3D_low, 
                              q3D_high)
#        print '2D XTAL: samplings = %5d, of which IN res range = %5d'%(N_samplings, 
#                                                                       N)            
        N_2D_exp_tot = N_2D_exp_tot + N
            
        # RECIPROCAL SPACE SPHERE
        N, N_samplings = getN(samplingSpace, 
                              qRod_limit_resSphere, 
                              q_2D, 
                              q3D_low, 
                              q3D_high)
#        print 'SPHERE:  samplings = %5d, of which IN res range = %5d'%(N_samplings, 
#                                                                       N)
        N_sphere_tot = N_sphere_tot + N
        
    ### Add ROD (0, 0) to N_sphere_tot ###
    q_2D = 0
    N, N_samplings = getN(samplingSpace, 
                          qLimit_3D, 
                          q_2D, 
                          q3D_low, 
                          q3D_high)
#    print 'SPHERE:  samplings = %5d, of which IN res range = %5d'%(N_samplings, 
#                                                                   N)
    N_sphere_tot = N_sphere_tot + N
    
    completeness = float(N_2D_exp_tot)/N_sphere_tot
                                                                
    return completeness
    

def calculateCompletenessPerBin(tiltAngle_deg,
                                cellSize,
                                wavelength,
                                d,
                                cutoff_2D,
                                overSampling):    
    
    # STANDARD PARAMETERS
    sinAlphas = numpy.linspace(-1, +1, 100, endpoint=True)
    hmax = 20
    kmax = 20
    
    tiltAngle = tiltAngle_deg * numpy.pi / 180 
    waveVector = 2 * numpy.pi / wavelength
    c_star = (2*numpy.pi)/(2*d)
    samplingSpace = c_star/overSampling
    qLimit_3D = 2*numpy.pi/binLimits[-1]
    
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(cellSize,
                                                                              hmax,
                                                                              kmax,
                                                                              cutoff_2D)
    completenesses = []
    nBins = len(binLimits)-1
    
    for i in range(0, nBins):
        print '\n\n\n\n********* BIN **********\n\n'
        
        dLimit_low  = binLimits[i]
        dLimit_high = binLimits[i+1]
        
        q3D_low  = 2*numpy.pi/dLimit_low
        q3D_high = 2*numpy.pi/dLimit_high
    
        print 'd [A]: %6.2f %6.2f    q [A-1]: %6.2f %6.2f\n\n'%(dLimit_low, 
                                                                dLimit_high,
                                                                q3D_low,
                                                                q3D_high)
                                                                
        completeness = getBinCompleteness(reciprocalLattice, 
                                          qLimit_3D, 
                                          sinAlphas, 
                                          tiltAngle,
                                          waveVector,
                                          samplingSpace,
                                          q3D_low, 
                                          q3D_high)
        
        completenesses.append(completeness)
        print 'Completeness:', completeness
        
        
    
       
    
    ### GLOBAL ###
    q3D_low  = 2*numpy.pi/binLimits[0]
    q3D_high = 2*numpy.pi/binLimits[-1]
    completeness = getBinCompleteness(reciprocalLattice, 
                                      qLimit_3D, 
                                      sinAlphas, 
                                      tiltAngle,
                                      waveVector,
                                      samplingSpace,
                                      q3D_low, 
                                      q3D_high)
                                      
    print '\n\nCOMPLETENESS PER BIN: ', completenesses    
    print 'GLOBAL COMPLETENESS: ', completeness


    
def calculateCompletenessPerTilt(cellSize,
                                 wavelength,
                                 d,
                                 cutoff_2D,
                                 overSampling,
                                 q3D_low, 
                                 q3D_high):   
                                     
    # STANDARD PARAMETERS
    tiltAngles_deg = [5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60.]
    sinAlphas = numpy.linspace(-1, +1, 100, endpoint=True)
    hmax = 20
    kmax = 20
        
    waveVector = 2 * numpy.pi / wavelength
    c_star = (2*numpy.pi)/(2*d)
    samplingSpace = c_star/overSampling
    qLimit_3D = q3D_high
    
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(cellSize,
                                                                              hmax,
                                                                              kmax,
                                                                              cutoff_2D)
    
    completenesses = []
    matplotlib.pyplot.figure(figsize=(10, 15))
    print 'Whole range (3D to 5.3A, consider only Bragg lines to 6A):'
    print 'TILT:   COMPLETENES:'
    for tiltAngle_deg in tiltAngles_deg:
        tiltAngle = tiltAngle_deg * numpy.pi / 180 
        completeness = getBinCompleteness(reciprocalLattice, 
                                          qLimit_3D, 
                                          sinAlphas, 
                                          tiltAngle,
                                          waveVector,
                                          samplingSpace,
                                          q3D_low, 
                                          q3D_high)
        completenesses.append(completeness)
        print tiltAngle_deg, completeness
            
    matplotlib.pyplot.figure()    
    matplotlib.pyplot.scatter(tiltAngles_deg, completenesses)
    axes = matplotlib.pyplot.gca()
    axes.tick_params(axis='x', labelsize=16)
    axes.tick_params(axis='y', labelsize=16)
    axes.set_xlabel(r'$\eta$ (degrees)', fontsize = 20) 
    axes.set_ylabel(r'Completeness',     fontsize = 20) 
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('./simulate_completeness.png', dpi=4*96)
        
    

if __name__ == "__main__":
    
    print "\n**** CALLING calculateCompletenessPerBin ****"
       
    tiltAngle_deg = 20.0 #deg
    cellSize = 62.45     #A
    wavelength = 1.48    #A
    d = 45               #A
    cutoff_2D = 6.0      #A
    overSampling = 2
    
    flag = 0
    if flag == 1:
        calculateCompletenessPerBin(tiltAngle_deg,
                                    cellSize,
                                    wavelength,
                                    d,
                                    cutoff_2D,
                                    overSampling)    
                                
                                
    q3D_low  = 2*numpy.pi/binLimits[0]
    q3D_high = 2*numpy.pi/binLimits[-1]
    
    flag = 1
    if flag == 1:
        calculateCompletenessPerTilt(cellSize,
                                     wavelength,
                                     d,
                                     cutoff_2D,
                                     overSampling,
                                     q3D_low,
                                     q3D_high)