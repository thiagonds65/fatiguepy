import numpy as np
import math
from . import prob_moment, Rainflow

class DK:
    def __new__(cls, *args, **kwargs):
        instance = super(DK, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        """
        The Acronym RR at the end of each method is equivalent to Rainflow Range Half Cycles
        and OR is equivalent to Ordinary Range Half Cycles
        """
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

    def PDFRR(self):
        
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

        # Dirlik method according to Matjaz
        #ps = ((G1/Q)*np.exp(-z/Q) + (G2*z/R**2)*np.exp(-z**2/(R**2))+G3*z*np.exp(-z**2/2))/(np.sqrt(self.m0))
        
        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += ps[i]*ds
        

        ps = ps/integ

        return ps

    def DamageRR(self):

        ps = self.PDFRR()
        ds = self.s[1] - self.s[0]
        DDK = 0
        for i in range(1,len(ps)):
            DDK += (self.EP*(self.C**(-1))*(self.s[i]**self.k)*ps[i]*ds)

        return DDK

    def LifesRR(self):
        TDKs = 1 / self.DamageRR()
        return TDKs

    def LifehRR(self):
        TDKh = self.LifesRR()/(3600)
        return TDKh

    def LifeRR(self):
        TDK = self.LifesRR()/self.xf
        return TDK
    
    def relative_errorRR(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            RR_value = self.LifeRR()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            RR_value = self.DamageRR()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(RR_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(RR_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
    
    def PDFOR(self):
        z = self.s / (np.sqrt(self.m0))
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        Xmin = (self.alpha2*(1+self.alpha2**2))/2
        G1 = (1/self.alpha2**2)*(Xm - Xmin)
        Q = 0.02+(2/self.alpha2)*(Xm - Xmin)
        G2 = 1 - (1/self.alpha2**2)*(Xm - Xmin)
        R = self.alpha2 + (1/self.alpha2)*(Xm - Xmin)
        psor = ((G1 / Q) * np.exp(-z / Q) + (G2 * z/ R ** 2) * np.exp(-z ** 2 / (2 * R ** 2)))/ (np.sqrt(self.m0))
        
        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += psor[i]*ds
        

        psor = psor/integ

        return psor

    def DamageOR(self):

        ps = self.PDFOR()
        ds = self.s[1] - self.s[0]
        DOR = 0
        for i in range(1,len(ps)):
            DOR += (self.EP*(self.C**(-1))*(self.s[i]**self.k)*ps[i]*ds)

        return DOR

    def LifesOR(self):
        TORs = 1 / self.DamageOR()
        return TORs

    def LifehOR(self):
        TORh = self.LifesOR()/(3600)
        return TORh

    def LifeOR(self):
        TOR = self.LifesOR()/self.xf
        return TOR
    
    def relative_errorOR(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            OR_value = self.LifeOR()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            OR_value = self.DamageOR()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(OR_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(OR_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
