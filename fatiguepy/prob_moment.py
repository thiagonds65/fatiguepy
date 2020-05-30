import numpy as np


class Probability_Moment:
    def __new__(cls, *args, **kwargs):
        instance = super(Probability_Moment, cls).__new__(cls)
        return instance

    def __init__(self, Y, f):
        self.Y = Y #Visto em Benasciutti Tovo sobre one-sided spectral density
        self.f = f
        self.df = f[1] - f[0]

    def moment0(self):

        m0 = 0

        # Por Dirlik, mi = (1/2pi) integral (w^i)*G(w)dw de -inf até inf

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m0 += self.Y[i]*self.df

        return m0

    def moment0dot75(self):
        m75 = 0

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m75 += np.abs(self.Y[i] * self.f[i]**0.75)*self.df

        return m75

    def moment1(self):
        m1 = 0

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m1 += np.abs(self.Y[i] * self.f[i])*self.df

        return m1

    def moment1dot5(self):
        m15 = 0

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m15 += np.abs(self.Y[i] * self.f[i]**1.5)*self.df

        return m15

    def moment2(self):
        m2 = 0


        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m2 += np.abs(self.Y[i] * self.f[i] ** 2)*self.df

        return m2

    def moment4(self):
        m4 = 0

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                m4 += np.abs(self.Y[i] * self.f[i] ** 4)*self.df

        return m4
    
    def E0(self):
        return np.sqrt(Probability_Moment(self.Y, self.f).moment2()/Probability_Moment(self.Y, self.f).moment0())
    
    def EP(self):
        return np.sqrt(Probability_Moment(self.Y, self.f).moment4()/Probability_Moment(self.Y, self.f).moment2())
    
    def alpha2(self):
        return Probability_Moment(self.Y, self.f).E0()/Probability_Moment(self.Y, self.f).EP()

