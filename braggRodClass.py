# -*- coding: utf-8 -*-
class braggRod:
    
    def __init__(self, hRod, kRod, qMin, qMax):
        self.hRod = hRod
        self.kRod = kRod
        self.qMin = qMin
        self.qMax = qMax      
        # braggRod objects objects are generated in merging


        
    def setExperimentalPoints(self, experimental_q, experimental_I):
        self.experimental_q = experimental_q
        self.experimental_I = experimental_I        
        
        
    def setModelPoints(self, model_q, model_I):
        self.model_q = model_q
        self.model_I = model_I
        
        
        
    def setModelCoefficients(self, coefficients):
        self.model_coefficients = coefficients
        
        
        
    def plotModel(self):
        import plotModel
        plotModel.plotModelFunction(self)


        
    def setExperimentalErrors(self, experimental_dI):
        self.experimental_dI = experimental_dI
        
        
        
    def setModelCoefficients_error(self, coefficients_error):
        self.model_coefficients_error = coefficients_error