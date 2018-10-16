# -*- coding: utf-8 -*-
import numpy

def sinc_function(x, a, N, delta, T, d):
    y = 0
    for i in range(-N, N+1):
        j = i+N
        if abs(x-(i*delta)) < 0.001: # Possible only if x is a single value
            y = y + a[j]
        else:
            y = (y + 
                 a[j]  * (numpy.sin(d*(x-i*delta)) / (d*(x-i*delta)) ) 
                       * numpy.exp(-0.5*(1/(4*numpy.pi**2))*T*(x-i*delta)**2))
    return y
    
def sigma_I_function(x, da, N, delta, T, d):
    var_I = 0
    for i in range(-N, N+1):
        j = i+N
        if abs(x-(i*delta)) < 0.001: # Possible only if x is a single value
            var_I = var_I + (da[j])**2
        else:
            var_I = (var_I + 
                      (da[j])**2  * 
                      (numpy.sin(d*(x-i*delta)) / (d*(x-i*delta)) )**2 *
                       numpy.exp(-(1/(4*numpy.pi**2))*T*(x-i*delta)**2)
                     )
    
    sigma_I = numpy.sqrt(var_I)
    
    return sigma_I