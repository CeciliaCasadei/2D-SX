# -*- coding: utf-8 -*-
import buildReciprocalLattice
import orbitClass

def makeOrbitsFunction(resolutionLimit):
    # BUILD RECIPROCAL LATTICE: h k qx qy resolution q
    cellSize = 62.45
    hmax = 100
    kmax = 100
    reciprocalLattice = buildReciprocalLattice.buildReciprocalLatticeFunction(cellSize, hmax, kmax, resolutionLimit)
    
    
    # ADD FLAG ITEM TO RECIPROCAL LATTICE SPOTS
    reciprocalLattice_flagged = []
    for reciprocalVector in reciprocalLattice:
        reciprocalVector_flagged = [reciprocalVector[0], reciprocalVector[1], 
                                    reciprocalVector[2], reciprocalVector[3], 
                                    reciprocalVector[4], reciprocalVector[5], 1]
        reciprocalLattice_flagged.append(reciprocalVector_flagged)
        
    # CHECK
    #    print '%5d %5d %10.4f %10.4f %12.2f %8.3f %d'%(reciprocalVector_flagged[0], reciprocalVector_flagged[1], 
    #                                                reciprocalVector_flagged[2], reciprocalVector_flagged[3], 
    #                                                reciprocalVector_flagged[4], reciprocalVector_flagged[5],
    #                                                reciprocalVector_flagged[6])
    #print len(reciprocalLattice_flagged)
    
    
    # BUILD ORBITS
    orbits = []
    for reflection_i in reciprocalLattice_flagged:
        h = reflection_i[0]
        k = reflection_i[1]
        flag = reflection_i[6]
        if flag == 1:
            reflection_i[6] = 0
            myOrbit = orbitClass.orbit(h, k)
            myOrbit.setResolution(reflection_i[4])
            for reflection_j in reciprocalLattice_flagged:
                if (reflection_j[0] == -h-k and reflection_j[1] == h) or (reflection_j[0] == k and reflection_j[1] == -h-k):
                    myOrbit.addOrbitIndices(reflection_j[0], reflection_j[1])
                    reflection_j[6] = 0   
            myOrbit.setOrbitLabel()
            orbits.append(myOrbit)
    
    # CHECK        
#    print len(orbits)
#    for orbit in orbits:
#        print orbit.orbitIndices, orbit.resolution, orbit.label
    
    return orbits