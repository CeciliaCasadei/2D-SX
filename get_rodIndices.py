# -*- coding: utf-8 -*-
import makeOrbits

def defineRodIndices(resolution_2D):
    orbits = makeOrbits.makeOrbitsFunction(resolution_2D)
    rodIndices = []
    for orbit in orbits:
        orbit_label = orbit.label
        if orbit_label[0] >= 0 and orbit_label[1] >= 0:
            rodIndices.append(orbit_label)  
    return rodIndices
