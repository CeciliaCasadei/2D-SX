# -*- coding: utf-8 -*-
import numpy

class orbit:
    
    def __init__(self, h, k):
        orbitIndices = []
        indicesToAdd = [h, k]
        orbitIndices.append(indicesToAdd)
        self.orbitIndices = orbitIndices
        

        
    def addOrbitIndices(self, h, k):
        indicesToAdd = [h, k]
        self.orbitIndices.append(indicesToAdd)
        
        
        
    def setResolution(self, resolution):
        self.resolution = resolution
        
        
        
    def setOrbitLabel(self):
        indices_1 = self.orbitIndices[0]
        indices_2 = self.orbitIndices[1]
        indices_3 = self.orbitIndices[2]
        h_1 = indices_1[0]
        k_1 = indices_1[1]
        h_2 = indices_2[0]
        k_2 = indices_2[1]
        h_3 = indices_3[0]
        k_3 = indices_3[1]
        
        # ASSIGN ORBIT LABEL:
        h_label = numpy.nan
        k_label = numpy.nan
        
        if (h_1 > 0 and k_1 > 0) or (h_1 < 0 and k_1 < 0) or (h_1 > 0 and k_1 == 0) or (h_1 < 0 and k_1 == 0):
            if not numpy.isnan(h_label):
                print "ERROR"
            if not numpy.isnan(k_label):
                print "ERROR"
                
            h_label = h_1
            k_label = k_1
            
        if (h_2 > 0 and k_2 > 0) or (h_2 < 0 and k_2 < 0) or (h_2 > 0 and k_2 == 0) or (h_2 < 0 and k_2 == 0):
            if not numpy.isnan(h_label):
                print "ERROR"
            if not numpy.isnan(k_label):
                print "ERROR"
                
            h_label = h_2
            k_label = k_2
            
        if (h_3 > 0 and k_3 > 0) or (h_3 < 0 and k_3 < 0) or (h_3 > 0 and k_3 == 0) or (h_3 < 0 and k_3 == 0):
            if not numpy.isnan(h_label):
                print "ERROR"
            if not numpy.isnan(k_label):
                print "ERROR"
                
            h_label = h_3
            k_label = k_3
                   
        self.label = [h_label, k_label]
    
        
        
