import numpy as np
import math
from . import prob_moment
class DK:
    def __new__(cls, *args, **kwargs):
        instance = super(DK, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.Y = Y
        self.f = f
        moments = prob_moment.Probability_Moment(self.Y, self.f)
        self.m0 = moments.moment0()
        self.m1 = moments.moment1()
        self.m2 = moments.moment2()
        self.m4 = moments.moment4()
        self.alpha2 = moments.alpha2()
        self.k = k
        self.C = C
        self.EP = moments.EP()
        self.xf = xf
        self.s = s
    
    def PDFOR(self):
        z = self.s / (2*np.sqrt(self.m0))
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        Xmin = (self.alpha2*(1+self.alpha2**2))/2
        G1 = (1/self.alpha2**2)*(Xm - Xmin)
        Q = 0.02+(2/self.alpha2)*(Xm - Xmin)
        G2 = 1 - (1/self.alpha2**2)*(Xm - Xmin)
        R = self.alpha2 + (1/self.alpha2)*(Xm - Xmin)
        psor = ((G1 / Q) * np.exp(-z / Q) + (G2 * z/ R ** 2) * np.exp(-z ** 2 / (2 * R ** 2)))/ (2 * np.sqrt(self.m0))
        
        return psor

    def PDF(self):
        z = self.s / (2*np.sqrt(self.m0))
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        G1 = 2 * (Xm - self.alpha2 ** 2) / (1 + self.alpha2 ** 2)
        R = (self.alpha2 - Xm - G1 ** 2) / (1 - self.alpha2 - G1 + G1 ** 2)
        G2 = (1 - self.alpha2 - G1 + G1 ** 2) / (1 - R)
        G3 = 1 - G1 - G2
        Q = 5 * (self.alpha2 - G3 - G2 * R) / (4 * G1)
        # # Abaixo contem a equacao do metodo Dirlik pela tese de mestrado de Ariduru
        p1 = (G1 / Q) * np.exp(-z / Q)
        p2 = (G2 / R ** 2) * np.exp(-z ** 2 / (R ** 2))
        p3 = G3 * z * np.exp(-z ** 2 / 2)
        pnum = p1 + p2 + p3
        pden = 2 * np.sqrt(self.m0)
        ps = pnum / pden
        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += ps[i]*ds

        ps = ps/integ

        # Abaixo contem a equacao do metodo Dirlik pelo artigo Matjaz
        # ps = ((G1/Q)*np.exp(-z/Q) + (G2*z/R**2)*np.exp(-z**2/(R**2))+G3*z*np.exp(-z**2/2))/(np.sqrt(self.m0))
        return ps

    def Damage(self):
        ps = self.PDF()
        ds = self.s[1] - self.s[0]
        DDK = 0
        for i in range(1,len(ps)):
            DDK += (self.EP*(self.C**(-1))*(self.s[i]**self.k)*ps[i]*ds)

        return DDK

    def Lifes(self):
        TDKs = 1 / self.Damage()
        return TDKs

    def Lifeh(self):
        TDKh = self.Lifes()/(3600)
        return TDKh

    def Life(self):
        TDK = self.Lifes()/self.xf
        return TDK