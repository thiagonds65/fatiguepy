import numpy as np
import math
from . import prob_moment, Rainflow

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
        if round(moments.alpha2(), 4) != 1.0:
            self.alpha2 = moments.alpha2()
        else:
            self.alpha2 = 0.99
        self.k = k
        self.C = C
        self.EP = moments.EP()
        self.xf = xf
        self.s = s

    def PDF(self):
        
        z = (self.s) / (np.sqrt(self.m0))
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        G1 = 2 * (Xm - self.alpha2 ** 2) / (1 + self.alpha2 ** 2)
        R = (self.alpha2 - Xm - G1 ** 2) / (1 - self.alpha2 - G1 + G1 ** 2)
        G2 = (1 - self.alpha2 - G1 + G1 ** 2) / (1 - R)
        G3 = 1 - G1 - G2
        Q = 5 * (self.alpha2 - G3 - G2 * R) / (4 * G1)

        # # Dirlik method according to Ariduru:
        p1 = (G1 / Q) * np.exp(-z / Q)
        p2 = (G2 * z / R ** 2) * np.exp(-z ** 2 / (2 * R ** 2))
        p3 = G3 * z * np.exp(-z ** 2 / 2)
        pnum = p1 + p2 + p3
        pden = np.sqrt(self.m0)
        ps = pnum / pden

        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += ps[i]*ds
        

        ps = ps/integ

        return ps
    
    def counting_cycles(self):
        ps = self.PDF()
        ds = self.s[1] - self.s[0]
        ns = ps*ds*self.EP*self.xf

        return ns
    
    def loading_spectrum(self):
        Cs = np.zeros(len(self.s))
        ns = self.counting_cycles()

        for i in range(len(self.s)):
            for j in range(i, len(self.s)):
                Cs[i] += ns[j]
        
        return Cs

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
    
    def relative_error(self, y, x, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            DK_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Life()
        elif type=="damage":
            DK_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = (DK_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = (DK_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
