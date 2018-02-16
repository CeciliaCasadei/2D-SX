# -*- coding: utf-8 -*-
import scipy.integrate
import numpy

# IMPLEMENTATION OF FRENCH-WILSON FUNCTION USING A STEP-WISE PRIOR AND GAUSSIAN p(I|J)

def Gauss_normalized(x, mu, sigma):
    y = (1/(sigma*numpy.sqrt(2*numpy.pi)))*numpy.exp(-0.5*((x-mu)**2)/(sigma**2))
    return y
    
def avg_intensity_integrand(x, mu, sigma):
    y = x * Gauss_normalized(x, mu, sigma)
    return y
    
def avg_amplitude_integrand(x, mu, sigma):
    y = numpy.sqrt(x) * Gauss_normalized(x, mu, sigma)
    return y
    
def var_intensity_integrand(x, mu, sigma):
    avg_intensity = scipy.integrate.quad(avg_intensity_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]
    y = (x-avg_intensity)**2 * Gauss_normalized(x, mu, sigma)
    return y
    
def var_amplitude_integrand(x, mu, sigma):
    avg_amplitude = scipy.integrate.quad(avg_amplitude_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]
    y = (numpy.sqrt(x)-avg_amplitude)**2 * Gauss_normalized(x, mu, sigma)
    return y
    
def calc_classic_std_amplitude(avg_intensity, std_intensity):
    if avg_intensity >= 0:
        return numpy.sqrt(avg_intensity+std_intensity)-numpy.sqrt(avg_intensity)
    else:
        return numpy.nan
        
def calc_avg_intensity(mu, sigma):
    return scipy.integrate.quad(avg_intensity_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]
    
def calc_avg_amplitude(mu, sigma):
    return scipy.integrate.quad(avg_amplitude_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]
    
def calc_var_intensity(mu, sigma):
    return scipy.integrate.quad(var_intensity_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]
    
def calc_var_amplitude(mu, sigma):
    return scipy.integrate.quad(var_amplitude_integrand, 0, mu+50*sigma, args=(mu, sigma))[0]

def testFunctions():    
    mu_sigmas = [(9., 2.), (10., 1.), (50., 6.), (0.5, 1.), (-0.4, 2.), (-5., 2.)]
    for mu_sigma in mu_sigmas:
        
        mu = mu_sigma[0]
        sigma = mu_sigma[1]
        
        normalization = scipy.integrate.quad(Gauss_normalized, mu-20*sigma, mu+20*sigma,  args=(mu, sigma))[0]
        avg_intensity = calc_avg_intensity(mu, sigma)
        avg_amplitude = calc_avg_amplitude(mu, sigma)
        var_intensity = calc_var_intensity(mu, sigma)    
        var_amplitude = calc_var_amplitude(mu, sigma)
        std_intensity = numpy.sqrt(var_intensity)
        std_amplitude = numpy.sqrt(var_amplitude)
        classic_std_amplitude = calc_classic_std_amplitude(mu, sigma)
        
        print '\n%.3f +- %.3f photons'%(mu, sigma)
        print 'NORMALIZATION: %.3f'%normalization
        print 'AVG INTENSITY: %.3f'%avg_intensity
        print 'AVG AMPLITUDE: %.3f'%avg_amplitude
        print 'VAR INTENSITY: %.3f'%var_intensity
        print 'STD INTENSITY: %.3f'%std_intensity
        print 'VAR AMPLITUDE: %.3f'%var_amplitude
        print 'STD AMPLITUDE: %.3f'%std_amplitude
        print 'CLASSIC STD AMPLITUDE: %.3f'%classic_std_amplitude
        
if __name__ == "__main__":
    testFunctions()
    