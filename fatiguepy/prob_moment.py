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

    def moment075(self):
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

    def moment15(self):
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
        return np.sqrt(self.moment2()/self.moment0())
    
    def EP(self):
        return np.sqrt(self.moment4()/self.moment2())
    
    def alpha1(self):
        return self.moment1()/np.sqrt(self.moment0()*self.moment2())

    def alpha2(self):
        return self.moment2()/np.sqrt(self.moment0()*self.moment4())
    
    def alpha075(self):
        return self.moment075()/np.sqrt(self.moment0()*self.moment15())

