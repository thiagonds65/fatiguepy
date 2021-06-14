import numpy as np


class Probability_Moment:
    def __new__(cls, *args, **kwargs):
        instance = super(Probability_Moment, cls).__new__(cls)
        return instance

    def __init__(self, Y, f):
        self.Y = Y #Visto em Benasciutti Tovo sobre one-sided spectral density
        self.f = f
        self.df = f[1] - f[0]

        '''
            If calculate m_i through integral of (G_yy(f)*(2*pi*f)^i)df, 
            the result of sqrt(m_2) is close to sigma_Xdot 
            and sqrt(m_4) is close to sigma_double dot, but the error is bigger than obtained
            without this (2*pi)
        '''

    def momentn(self, n):

        mn = 0

        # Por Dirlik, mi = (1/2pi) integral (w^i)*G(w)dw de -inf até inf

        for i in range(0, len(self.f)):
            if self.f[i] >= 0:
                mn += np.abs(self.Y[i]*(2*np.pi*self.f[i])**n)*self.df

        return mn
    
    def E0(self):
        return np.sqrt(self.momentn(2)/self.momentn(0))/(2*np.pi)
    
    def EP(self):
        return np.sqrt(self.momentn(4)/self.momentn(2))/(2*np.pi)
    
    def alphan(self, n):
        return self.momentn(n)/np.sqrt(self.momentn(0)*self.momentn(2*n))


