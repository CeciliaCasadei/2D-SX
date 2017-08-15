# -*- coding: utf-8 -*-
class partialSum:
    
    def __init__(self, h_label, k_label, module_bottomBound, module_topBound, module_leftBound, module_rightBound, nTerms, partialSum):
        self.h_label = h_label
        self.k_label = k_label
        self.module_bottomBound = module_bottomBound
        self.module_topBound = module_topBound
        self.module_leftBound = module_leftBound
        self.module_rightBound = module_rightBound
        self.nTerms = nTerms
        self.partialSum = partialSum # NOT BG SUBTRACTED, NOT NORMALIZED, N OF COUNTS