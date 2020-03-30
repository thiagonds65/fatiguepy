import math
import numpy as np

class NB:
    def __new__(cls, *args, **kwargs):
        instance = super(NB, cls).__new__(cls)
        return instance

    def __init__(self, k, E0, C, m0, xf):
        self.E0 = E0
        self.k = k
        self.C = C
        self.m0 = m0
        self.xf = xf

    def Damage(self):
        EGFkbar2 = math.gamma(1 + self.k / 2)  # Euler Gamma Function

        DNB = self.E0 * (self.C ** (-1)) * (np.sqrt(2 * self.m0)) ** self.k * EGFkbar2 * self.xf
        return DNB
    
    def Life(self):
        TNB = 1/NB(self.k, self.E0, self.C, self.m0, self.xf).Damage()
        return TNB
    
    def Lifeh(self):
        TNBh = self.xf * NB(self.k, self.E0, self.C, self.m0, self.xf).Life() / 3600
        return TNBh
