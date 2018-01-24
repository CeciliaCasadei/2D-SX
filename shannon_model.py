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
    
