# -*- coding: utf-8 -*-
import numpy
a = 62.45
b = 62.45
c = 100.9
directVectors_cartesian = numpy.matrix([[a, b*numpy.cos(2*numpy.pi/3), 0], [0, b*numpy.sin(2*numpy.pi/3), 0], [0, 0, c]])

p3_cartesian = numpy.matrix([[numpy.cos(2*numpy.pi/3), -numpy.sin(2*numpy.pi/3), 0],[numpy.sin(2*numpy.pi/3), numpy.cos(2*numpy.pi/3), 0], [0, 0, 1]]) 
p3_crystallographic = directVectors_cartesian.I * p3_cartesian * directVectors_cartesian
print p3_crystallographic

p3p3_cartesian = numpy.matrix([[numpy.cos(4*numpy.pi/3), -numpy.sin(4*numpy.pi/3), 0],[numpy.sin(4*numpy.pi/3), numpy.cos(4*numpy.pi/3), 0], [0, 0, 1]]) 
p3p3_crystallographic = directVectors_cartesian.I * p3p3_cartesian * directVectors_cartesian
print p3p3_crystallographic

cartesianToCrystallographic = directVectors_cartesian.I 
print cartesianToCrystallographic

reciprocalVectors_cartesian = 2 * numpy.pi * directVectors_cartesian.I.T
print reciprocalVectors_cartesian
