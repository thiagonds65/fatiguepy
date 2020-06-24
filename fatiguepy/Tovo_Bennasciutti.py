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
        moment = prob_moment.Probability_Moment(self.Y, self.f)
        self.alpha1 = prob_moment.Probability_Moment(self.Y, self.f).alpha1()
        self.alpha2 = prob_moment.Probability_Moment(self.Y, self.f).alpha2()
        self.E0 = prob_moment.Probability_Moment(self.Y, self.f).E0()
        self.DNB = Narrow_Band.NB(self.k, self.C, self.Y, self.f, self.xf, self.s).Damage()

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
    
    def relative_error(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            TB_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            TB_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(TB_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(TB_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
