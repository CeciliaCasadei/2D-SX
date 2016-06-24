# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 14:42:24 2016
@author: casadei_c
MEMBER OF Lattice CLASS
PRODUCE TEXT FILE CONTAINING THE REFINED PREDICTED PATTERN
"""

def logRefinedPatternInfoFunction(self, folderName):
    fSummary = open('%s/r%s_Image%s_Lattice%s.txt'%(folderName, self.runNumber, self.imageNumber, self.latticeNumberInImage), 'wb')  
    fSummary.write('Run: %s \nImage number: %s \nLattice number: %s'%(self.runNumber, self.imageNumber, self.latticeNumberInImage)) 
    fSummary.write('\nImage filename: %s \nTilt angle: %s \nOriginal in-plane rotation: %.3f \nRefined in-plane rotation: %.3f'%(self.fileName, 
                                                                                                                                 self.tiltAngle, 
                                                                                                                                 self.inPlaneRotation, 
                                                                                                                                 self.refinedInPlaneOrientation))
    fSummary.write('\nOriginal cell size: %.2f \nRefined cell size: %.2f'%(self.referenceCellSize, self.refinedCellSize))
    fSummary.write('\nRefined predicted pattern:\n')
    fSummary.write('  h    k        q_x        q_y      d_min        q   azimuth    rotated azimuth    azimuth on detector    diffraction angle    radius on detector      qRod    LP\n')
    for predictedPatternLine in self.refinedPredictedPattern:
        fSummary.write('%3d%5d%11.4f%11.4f%11.2f%9.3f%10.3f%19.4f%23.4f%21.4f%22.4f%10.5f%6.2f\n'
                        %(predictedPatternLine[0], predictedPatternLine[1], 
                          predictedPatternLine[2], predictedPatternLine[3], 
                          predictedPatternLine[4], predictedPatternLine[5], 
                          predictedPatternLine[6], predictedPatternLine[7], 
                          predictedPatternLine[8], predictedPatternLine[9], 
                          predictedPatternLine[10], predictedPatternLine[11], 
                          predictedPatternLine[12]))
    fSummary.close() 