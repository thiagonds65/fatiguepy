import numpy as np
from . import prob_moment, Narrow_Band, Rainflow

class AL:
    def __new__(cls, *args, **kwargs):
        instance = super(AL, cls).__new__(cls)
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
        self.alpha075 = moments.alphan(0.75)
        self.alpha1 = moments.alphan(1)
        self.alpha2 = moments.alphan(2)
        self.EP = moments.EP()
        NB = Narrow_Band.NB(self.k, self.C, self.Y, self.f, self.xf, self.s)
        self.DNB = NB.Damage()
        self.pNB = NB.PDF()

    def PDF(self):

        lambdaAL = self.alpha075**2
        pAL = lambdaAL*self.pNB

        return pAL
    
    def counting_cycles(self):
        pAL = self.PDF()
        ds = self.s[1] - self.s[0]
        nAL = pAL*ds*self.EP*self.xf

        return nAL
    
    def loading_spectrum(self):
        CAL = np.zeros(len(self.s))
        nAL = self.counting_cycles()

        for i in range(len(self.s)):
            for j in range(i, len(self.s)):
                CAL[i] += nAL[j]
        
        return CAL

    def Damage(self):

        DAL = (self.alpha075**2) * self.DNB

        return DAL

    def Lifes(self):  # ciclos
        TALs = 1 / self.Damage()
        return TALs

    def Lifeh(self):  # horas
        TALh = self.xf * self.Life() / (60 * 60)
        return TALh
    
    def Life(self):
        TAL = self.Lifes()/self.xf
        return TAL
    
    def relative_error(self, y, x, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            AL_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Life()
        elif type=="damage":
            AL_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = (AL_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = (AL_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
