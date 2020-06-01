import numpy as np
import math
from . import prob_moment, Rainflow

class ZB:
    def __new__(cls, *args, **kwargs):
        instance = super(ZB, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s, zb):
        self.Y = Y
        self.f = f
        moments = prob_moment.Probability_Moment(self.Y, self.f)
        self.m0 = moments.moment0()
        self.m1 = moments.moment1()
        self.m2 = moments.moment2()
        self.m4 = moments.moment4()
        self.m75 = moments.moment0dot75()
        self.m15 = moments.moment1dot5()
        self.alpha2 = moments.alpha2()
        self.k = k
        self.C = C
        self.EP = moments.EP()
        self.xf = xf
        self.zb = zb
        self.s = s

    def beta(self):
        if self.alpha2<0.9:
            beta = 1.1
        else:
            beta = 1.1+9*(self.alpha2-0.9)
        return beta


    def PDF(self):
        z = self.s / (2*np.sqrt(self.m0))
        beta = ZB(self.k, self.C, self.Y, self.f, self.xf, self.s, self.zb).beta()
        EGFbeta = math.gamma(1 + 1/beta)
        EGF3beta = math.gamma(1+3/beta)

        if self.zb==1:
            alpha = 8 - 7*self.alpha2 #for ZB1
        elif self.zb==2:
            
            alpha75 = self.m75/(np.sqrt(self.m0*self.m15))
            print(f"alpha75 = {alpha75}")
            if alpha75 < 0.5:
                roZB = 0.28
            else:
                roZB = -0.4154 + 1.392*self.alpha2
            d = 0
            for i in range(10):
                d = d - (EGF3beta*(1-self.alpha2**2)*d**3+3*EGFbeta*(roZB*self.alpha2-1)*d+3*np.sqrt(np.pi/2)*self.alpha2*(1-roZB))\
                    /(3*EGF3beta*(1-self.alpha2**2)*d**2+3*EGFbeta*(roZB*self.alpha2-1))
                alpha = d**(-beta)

        w = (1 - self.alpha2)/(1 - math.sqrt(2/math.pi)*EGFbeta*alpha**(-1/beta))

        ps = (w*alpha*beta*z**(beta - 1)*np.exp(-alpha*z**beta) + (1-w)*z*np.exp((-z**2)/2)/(2*np.sqrt(self.m0))) 
        
        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += ps[i]*ds
        ps = ps/integ

        return ps

    def Damage(self):
        ps = self.PDF()
        ds = self.s[1] - self.s[0]
        DZB = 0
        for i in range(1,len(ps)):
            DZB += abs(self.EP*(self.C**(-1))*(self.s[i]**self.k)*ps[i]*ds)

        return DZB

    def Lifes(self):
        TZB = 1 / self.Damage()
        return TZB

    def Lifeh(self):
        TZBh = self.Lifes()/(3600)
        return TZBh

    def Life(self):
        TZB = self.Lifes()/self.xf
        return TZB

    def relative_error(self, y, method="Rainflow", Dexperimental=None):
        DZB = self.Damage()
        if(method == "Rainflow"):
            DRF = Rainflow.rainflowD(self.C, self.k, y, self.xf).DRF()
            err = abs(DZB - DRF)/DRF
        elif(method == "Experimental" and Dexperimental != None):
            DEX = Dexperimental
            err = abs(DZB - DEX)/DEX
        elif(method == "Experimental" and Dexperimental == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
