import math
import numpy as np
from . import prob_moment, Rainflow
from scipy.stats import norm

class RC:
    def __new__(cls, *args, **kwargs):
        instance = super(RC, cls).__new__(cls)
        return instance

    def __init__(self, k, C, Y, f, xf, s):
        self.k = k
        self.C = C
        self.xf = xf
        self.Y = Y
        self.f = f
        self.s = s
        moments = prob_moment.Probability_Moment(self.Y, self.f)
        self.m0 = moments.moment0()
        self.m1 = moments.moment1()
        self.m2 = moments.moment2()
        self.E0 = moments.E0()
        self.alpha2 = moments.alpha2()
        self.EP = moments.EP()

    def PDF(self):
        """
        A method to obtain Probability Density Function based in Rice's PDF Process.

        Parameters
        ----------

        """
        
        
        if round(self.alpha2, 4) != 1.0:
            z = self.s/np.sqrt(self.m0)
            ratio = (np.sqrt(1-self.alpha2**2)/np.sqrt(2*np.pi*self.m0))
            exp = np.exp(-(z**2)/(2*(1-self.alpha2**2)))
            
            ratio2 = (self.alpha2*self.s/self.m0)
            exp2 = np.exp(-(z**2)/2)

            #Error Function
            '''
            Considering phi the standard normal distribution and erf the error function,
            the relation between this two functions is given by:
            phi(x) = (1/2)*(1+erf(x/sqrt(2)))
            or
            erf(x) = 2*phi(x*sqrt(2))-1
            '''
            x = self.alpha2*z/np.sqrt(1-self.alpha2**2)
            for i in range(len(self.s)):
                phi = (1/2)*(math.erf(x[i]/np.sqrt(2))+1)

            # # Abaixo contem a função de distribuição normal pelo artigo de Carpinteri
            pRC = ratio*exp + ratio2*exp2*phi
        else:
            z = self.s/np.sqrt(self.m0)
            ratio = (self.s/self.m0)
            exp = np.exp(-(z**2)/2)
        
            pRC = ratio*exp
        
        return pRC
    
    def counting_cycles(self):
        pRC = self.PDF()
        ds = self.s[1] - self.s[0]
        nRC = pRC*ds*self.EP*self.xf

        return nRC
    
    def loading_spectrum(self):
        CRC = np.zeros(len(self.s))
        nRC = self.counting_cycles()

        for i in range(len(self.s)):
            for j in range(i, len(self.s)):
                CRC[i] += nRC[j]
        
        return CRC

    def Damage(self):

        pRC = self.PDF()
        ds = self.s[1] - self.s[0]
        DRC = 0
        for i in range(1,len(pRC)):
            DRC += self.EP*(self.C**(-1))*(self.s[i]**self.k)*pRC[i]*ds

        return DRC
    
    def Lifes(self):
        TRCs = 1/self.Damage()
        return TRCs
    
    def Lifeh(self):
        TRCh = self.Lifes()/3600
        return TRCh

    def Life(self):
        TRC = self.Lifes()/self.xf
        return TRC
    
    def relative_error(self, y, x, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            RC_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Life()
        elif type=="damage":
            RC_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, x).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = (RC_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = (RC_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
