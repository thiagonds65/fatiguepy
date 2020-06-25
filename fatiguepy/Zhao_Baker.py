import numpy as np
import math
from . import prob_moment, Rainflow

class ZB:
    def __new__(cls, *args, **kwargs):
        instance = super(ZB, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.Y = Y
        self.f = f
        moments = prob_moment.Probability_Moment(self.Y, self.f)
        self.m0 = moments.moment0()
        self.m1 = moments.moment1()
        self.m2 = moments.moment2()
        self.m4 = moments.moment4()
        self.m75 = moments.moment075()
        self.m15 = moments.moment15()
        self.alpha2 = moments.alpha2()
        self.k = k
        self.C = C
        self.EP = moments.EP()
        self.xf = xf
        self.s = s

    def beta(self):
        if self.alpha2<0.9:
            beta = 1.1
        else:
            beta = 1.1+9*(self.alpha2-0.9)
        return beta


    def PDF(self):
        z = self.s / (np.sqrt(self.m0))
        beta = self.beta()
        EGFbeta = math.gamma(1 + 1/beta)
        EGF3beta = math.gamma(1+3/beta)

        if self.k>3:
            alpha = 8 - 7*self.alpha2 #for ZB1
        elif self.k<=3:
            
            alpha75 = self.m75/(np.sqrt(self.m0*self.m15))

            if alpha75 < 0.5:
                roZB = 0.28
            else:
                roZB = -0.4154 + 1.392*alpha75
            d = 10
            for i in range(10):
                d = d - (EGF3beta*(1-self.alpha2**2)*d**3+3*EGFbeta*(roZB*self.alpha2-1)*d+3*np.sqrt(np.pi/2)*self.alpha2*(1-roZB))\
                    /(3*EGF3beta*(1-self.alpha2**2)*d**2+3*EGFbeta*(roZB*self.alpha2-1))
                alpha = d**(-beta)
        

        w = (1 - self.alpha2)/(1 - math.sqrt(2/math.pi)*EGFbeta*alpha**(-1/beta))

        ps = (w*alpha*beta*z**(beta - 1)*np.exp(-alpha*z**beta) + (1-w)*z*np.exp((-z**2)/2))/(np.sqrt(self.m0)) 
        
        '''
        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += ps[i]*ds
        
        print(integ)  

        denom = self.EP*self.xf # First hypothesis
        ps = ps/integ
        '''

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

    def relative_error(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            ZB_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            ZB_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(ZB_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(ZB_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
