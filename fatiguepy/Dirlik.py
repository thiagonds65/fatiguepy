import numpy as np
import math
class DK:
    def __new__(cls, *args, **kwargs):
        instance = super(DK, cls).__new__(cls)
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

    def PDF(self):
        z = self.s / (2*np.sqrt(self.m0))
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        G1 = 2 * (Xm - self.alpha2 ** 2) / (1 + self.alpha2 ** 2)
        R = (self.alpha2 - Xm - G1 ** 2) / (1 - self.alpha2 - G1 + G1 ** 2)
        G2 = (1 - self.alpha2 - G1 + G1 ** 2) / (1 - R)
        G3 = 1 - G1 - G2
        Q = 5 * (self.alpha2 - G1 - G2 * R) / (4 * G1)
        # # Abaixo contem a equacao do metodo Dirlik pela tese de mestrado de Ariduru
        ps = ((G1 / Q) * np.exp(-z / Q) + (G2 / R ** 2) * np.exp(-z ** 2 / (R ** 2)) \
        + G3 * z * np.exp(-z ** 2 / 2)) / (2 * np.sqrt(self.m0))
        # Abaixo contem a equacao do metodo Dirlik pelo artigo Matjaz
        # ps = ((G1/Q)*np.exp(-z/Q) + (G2*z/R**2)*np.exp(-z**2/(R**2))+G3*z*np.exp(-z**2/2))/(np.sqrt(self.m0))
        return ps

    def Damage(self):
        EGFk = math.gamma(1 + self.k)
        Xm = (self.m1 / self.m0) * np.sqrt(self.m2 / self.m4)
        G1 = 2 * (Xm - self.alpha2 ** 2) / (1 + self.alpha2 ** 2)
        R = (self.alpha2 - Xm - G1 ** 2) / (1 - self.alpha2 - G1 + G1 ** 2)
        G2 = (1 - self.alpha2 - G1 + G1 ** 2) / (1 - R)
        G3 = 1 - G1 - G2
        Q = 5 * (self.alpha2 - G1 - G2 * R) / (4 * G1)
        DDK = self.C ** (-1) * self.EP * self.m0 ** (self.k / 2) * (
                G1 * (Q ** self.k) * EGFk + (np.sqrt(2)) **
                self.k * math.gamma(self.k / 2 + 1) * (G2 * np.abs(R) ** self.k + G3))
        return DDK

    def Life(self):
        TDK = 1 / DK(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).Damage()
        return TDK

    def Lifeh(self):
        TDKh = self.xf * \
               DK(self.m0, self.m1, self.m2, self.m4, self.alpha2, self.s, self.k, self.C, self.EP, self.xf).Life() \
               / (60 * 60)
        return TDKh
