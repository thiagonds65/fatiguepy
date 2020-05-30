import math
import numpy as np
from . import prob_moment
from scipy.stats import norm

class NB:
    def __new__(cls, *args, **kwargs):
        instance = super(NB, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.k = k
        self.C = C
        self.xf = xf
        self.Y = Y
        self.f = f
        self.s = s
        self.m0 = prob_moment.Probability_Moment(self.Y, self.f).moment0()
        self.E0 = prob_moment.Probability_Moment(self.Y, self.f).E0()
        self.alpha2 = prob_moment.Probability_Moment(self.Y, self.f).alpha2()
        self.EP = prob_moment.Probability_Moment(self.Y, self.f).EP()

    def PDPeaks(self):
        ratio = (np.sqrt(1-self.alpha2**2)/np.sqrt(2*np.pi*self.m0))
        exp = np.exp(-(self.s**2)/(2*self.m0*(1-self.alpha2**2)))
        ratio2 = (self.alpha2*self.s/self.m0)
        exp2 = np.exp(-(self.s**2)/(2*self.m0))

        #Error Function
        '''
        Considering phi the standard normal distribution and erf the error function,
        the relation between this two functions is given by:
        phi(x) = (1/2)*(1+erf(x/sqrt(2)))
        or
        erf(x) = 2*phi(x*sqrt(2))-1
        '''

        x = (self.alpha2*self.s)/np.sqrt(self.m0*(1-self.alpha2**2))
        for i in range(len(self.s)):
            phi = (1/2)*(math.erf(x[i]/np.sqrt(2))+1)
        
        z = self.s / (np.sqrt(self.m0))
        # # Abaixo contem a função de distribuição normal pelo artigo de Carpinteri
        #pp = ratio*exp + ratio2*exp2*phi
        # Visto no artigo Dirlik (página 63 do pdf):
        pp = (self.s/self.m0)*np.exp(-self.s**2/(2*self.m0))

        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += pp[i]*ds
        pp = pp/integ
        
        return pp

    def Damage(self):
        '''
        O dano calculado pela aproximação da banda estreita pode ser expresso 
        pela seguinte expressão (MRSNIK et. al, 2016)
        EGFkbar2 = math.gamma(1 + self.k / 2)  # Euler Gamma Function

        DNB = self.E0 * (self.C ** (-1)) * (np.sqrt(2 * self.m0)) ** self.k * EGFkbar2 * self.xf
        Mas também pode ser expresso através de uma equação empírica dependente da PDF, 
        como visto abaixo:
        '''
        pp = self.PDPeaks()
        ds = self.s[1] - self.s[0]
        DNB = 0
        for i in range(1,len(pp)):
            DNB += self.EP*(self.C**(-1))*(self.s[i]**self.k)*pp[i]*ds

        return DNB
    
    def Lifes(self):
        TNBs = 1/self.Damage()
        return TNBs
    
    def Lifeh(self):
        TNBh = self.Lifes()/3600
        return TNBh

    def Life(self):
        TNB = self.Lifes()/self.xf
        return TNB
