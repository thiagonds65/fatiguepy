import numpy as np
from . import prob_moment, Narrow_Band, Rainflow

class WL:
    def __new__(cls, *args, **kwargs):
        instance = super(WL, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.k = k
        self.C = C
        self.Y = Y
        self.f = f
        self.xf = xf
        self.s = s
        moment = prob_moment.Probability_Moment(self.Y, self.f)
        self.alpha2 = prob_moment.Probability_Moment(self.Y, self.f).alpha2()
        self.E0 = prob_moment.Probability_Moment(self.Y, self.f).E0()
        self.DNB = Narrow_Band.NB(self.k, self.C, self.Y, self.f, self.xf, self.s).Damage()

    def Damage(self):
        epsylon = np.sqrt(np.abs(1 - self.alpha2 ** 2))

        ak = 0.926 - 0.033 * self.k
        bk = 1.587 * self.k - 2.323

        roWL = ak + (1 - ak) * (1 - epsylon) ** bk

        DWL = roWL * self.DNB

        return DWL

    def Lifes(self):  # ciclos
        TWLs = 1 / self.Damage()
        return TWLs

    def Lifeh(self):  # horas
        TWLh = self.xf * self.Life() / (60 * 60)
        return TWLh
    
    def Life(self):
        TWL = self.Lifes()/self.xf
        return TWL
    
    def relative_error(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            WL_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            WL_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(WL_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(WL_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
