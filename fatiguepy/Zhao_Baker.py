import numpy as np
import math
class ZB:
    def __new__(cls, *args, **kwargs):
        instance = super(ZB, cls).__new__(cls)
        return instance

    def __init__(self, m0, m1, m2, m4, alpha2, s, k, C, EP, xf):
        self.m0 = m0
        self.m1 = m1
        self.m2 = m2
        self.m4 = m4
        self.alpha2 = alpha2
        self.s = s
        self.k = k
        self.C = C
        self.EP = EP
        self.xf = xf

    def beta(self):
        if self.alpha2<0.9:
            beta = 1.1
        else:
            beta = 1.1+9*(self.alpha2-0.9)
        return beta


    def PDF(self):
        z = self.s / (2*np.sqrt(self.m0))
        beta = ZB(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).beta()
        alpha = 8 - 7*self.alpha2
        EGFbeta = math.gamma(1 + 1/beta)
        w = (1 - self.alpha2)/(1 - math.sqrt(2/math.pi)*EGFbeta*alpha**(-1/beta))

        ps = w*alpha*beta*z**(beta - 1)*np.exp(-alpha*z**beta) + (1-w)*z*np.exp((-z**2)/2) 
        # Abaixo contem a equacao do metodo Dirlik pelo artigo Matjaz
        # ps = ((G1/Q)*np.exp(-z/Q) + (G2*z/R**2)*np.exp(-z**2/(R**2))+G3*z*np.exp(-z**2/2))/(np.sqrt(self.m0))
        return ps

    def Damage(self):
        beta = ZB(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).beta()
        alpha = 8 - 7*self.alpha2
        EGFbeta = math.gamma(1 + 1/beta)
        EGFk = math.gamma(1 + self.k/beta)
        EGFk2 = math.gamma(1 + self.k/2)
        w = (1 - self.alpha2)/(1 - math.sqrt(2/math.pi)*EGFbeta*alpha**(-1/beta))

        DZB = (self.EP/self.C)*self.m0**(self.k/2)*(w*alpha**(-self.k/beta)*EGFk+(1-w)*2**(self.k/2)*EGFk2)
        return DZB

    def Life(self):
        TZB = 1 / ZB(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).Damage()
        return TZB

    def Lifeh(self):
        TZBh = self.xf * \
               ZB(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).Life() \
               / (60 * 60)
        return TZBh