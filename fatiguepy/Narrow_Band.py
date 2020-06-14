import math
import numpy as np
from . import prob_moment, Rainflow
from scipy.stats import norm

class NB:
    def __new__(cls, *args, **kwargs):
        instance = super(NB, cls).__new__(cls)
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

    def PDPeaks(self):
        """
        A method to obtain Probability Density Function for Narrow Band Process.

        Parameters
        ----------

        """
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
        
        z = self.s / (np.sqrt(self.m0))
        # # Abaixo contem a função de distribuição normal pelo artigo de Carpinteri
        pp = ratio*exp + ratio2*exp2*phi

        integ = 0
        ds = self.s[1] - self.s[0]
        for i in range(len(self.s)):
            integ += pp[i]*ds
        pp = pp/integ
        
        return pp

    def Damage(self):
        '''
        O dano calculado pela aproximação da banda estreita pode ser expresso 
        pela seguinte expressão (MRSNIK et. al, 2016)
        EGFkbar2 = math.gamma(1 + self.k / 2)  # Euler Gamma Function

        DNB = self.E0 * (self.C ** (-1)) * (np.sqrt(2 * self.m0)) ** self.k * EGFkbar2 * self.xf
        Mas também pode ser expresso através de uma equação empírica dependente da PDF, 
        como visto abaixo:
        '''
        pp = self.PDPeaks()
        ds = self.s[1] - self.s[0]
        DNB = 0
        for i in range(1,len(pp)):
            DNB += self.EP*(self.C**(-1))*(self.s[i]**self.k)*pp[i]*ds
        
        #DNB = self.E0*self.C**(-1)*(2*self.m0)**(self.k/2)*math.gamma(1+self.k/2)

        return DNB
    
    def Lifes(self):
        TNBs = 1/self.Damage()
        return TNBs
    
    def Lifeh(self):
        TNBh = self.Lifes()/3600
        return TNBh

    def Life(self):
        TNB = self.Lifes()/self.xf
        return TNB
    
    def relative_error(self, y, method="Rainflow", experimental_value=None, type='cycles'):
        if type=="cycles":
            NB_value = self.Life()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Life()
        elif type=="damage":
            NB_value = self.Damage()
            RF_value = Rainflow.rainflowD(self.C, self.k, y, self.xf).Damage()
        elif type!="cycles" and type!="damage":
            raise UnboundLocalError("Invalid type. Try 'cycles' or 'damage'")
        
        if(method == "Rainflow"):
            err = abs(NB_value - RF_value)/RF_value
        elif(method == "Experimental" and experimental_value != None):
            EX_value = experimental_value
            err = abs(NB_value - EX_value)/EX_value
        elif(method == "Experimental" and experimental_value == None):
            raise UnboundLocalError("Dexperimental must be different from None for method 'Experimental'")
        elif(method != "Experimental" and method != "Rainflow"):
            raise UnboundLocalError("Invalid Method. Try method='Rainflow' or method='Experimental'")

        return err
