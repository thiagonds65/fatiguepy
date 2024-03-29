import numpy as np
from . import prob_moment, Narrow_Band, Rainflow

class TB:
    def __new__(cls, *args, **kwargs):
        instance = super(TB, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.k = k
        self.C = C
        self.Y = Y
        self.f = f
        self.xf = xf
        self.s = s
        moments = prob_moment.Probability_Moment(self.Y, self.f)
        self.m0 = moments.momentn(0)
        self.m1 = moments.momentn(1)
        self.m2 = moments.momentn(2)
        self.E0 = moments.E0()
        self.EP = moments.EP()
        self.alpha1 = moments.alphan(1)
        self.alpha2 = moments.alphan(2)
        NB = Narrow_Band.NB(self.k, self.C, self.Y, self.f, self.xf, self.s)
        self.DNB = NB.Damage()
        self.pNB = NB.PDF()
    
    def PDF(self):
        parameter = (self.alpha1 - self.alpha2)/(1 - self.alpha1)
        b = min([1, parameter])

        lambdaTB = (b + (1 - b) * self.alpha2**(self.k - 1))*self.alpha2
        pTB = lambdaTB*self.pNB

        # pTB = (self.s/(self.m0*self.alpha2**2))*np.exp(-self.s**2/(2*self.m0*self.alpha2**2))

        return pTB
    
    def counting_cycles(self):
        pTB = self.PDF()
        ds = self.s[1] - self.s[0]
        nTB = pTB*ds*self.EP*self.xf

        return nTB
    
    def loading_spectrum(self):
        CTB = np.zeros(len(self.s))
        nTB = self.counting_cycles()

        for i in range(len(self.s)):
            for j in range(i, len(self.s)):
                CTB[i] += nTB[j]
        
        return CTB

    def Damage(self):
        parameter = (self.alpha1 - self.alpha2)/(1 - self.alpha1)
        b = min([1, parameter])
        '''
        diff = self.alpha1 - self.alpha2
        eq = (1 + self.alpha1*self.alpha2 - (self.alpha1 + self.alpha2))*np.exp(2.11*self.alpha2)
        denom = (self.alpha2 - 1)**2
        b = diff*(1.112*eq + diff)/denom
        '''
        DTB = (b + (1 - b) * self.alpha2**(self.k - 1)) * self.alpha2 * self.DNB

        return DTB

    def Lifes(self):  # ciclos
        TTBs = 1 / self.Damage()
        return TTBs

    def Lifeh(self):  # horas
        TTBh = self.xf * self.Life() / (60 * 60)
        return TTBh
    
    def Life(self):
        TTB = self.Lifes()/self.xf
        return TTB
    
    def relative_error(self, y, x, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            TB_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Life()
        elif type=="damage":
            TB_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = (TB_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = (TB_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
