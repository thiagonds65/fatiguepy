import numpy as np
import math
class prob_dist:
    def __new__(cls, *args, **kwargs):
        instance = super(prob_dist, cls).__new__(cls)
        return instance

    def __init__(self, m0, m1, m2, m4, alpha2, s, k, C, EP, xf):
        self.m0 = m0
        self.m1 = m1
        self.m2 = m2
        self.m4 = m4
        self.alpha2 = alpha2
        self.s = s
        self.k = k
        self.C = C
        self.EP = EP
        self.xf = xf

    def PDPeaks(self):
        ratio = (np.sqrt(1-self.alpha2**2)/np.sqrt(2*np.pi*self.m0))
        exp = np.exp(-(self.s**2)/(2*self.m0*(1-self.alpha2**2)))
        ratio2 = (self.alpha2*self.s/self.m0)
        exp2 = np.exp(-(self.s**2)/(2*self.m0))
        #Error Function
        x = (self.alpha2*self.s)/np.sqrt(self.m0*(1-self.alpha2**2))
        for i in range(len(self.s)):
            phi = (1/2)*(math.erf(x[i]/np.sqrt(2))+1)
        
        z = self.s / (np.sqrt(self.m0))
        # # Abaixo contem a função de distribuição normal pelo artigo de Carpinteri
        pp = ratio*exp + ratio2*exp2*phi
        return pp

