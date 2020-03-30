import numpy as np


class Probability_Moment:
    def __new__(cls, *args, **kwargs):
        instance = super(Probability_Moment, cls).__new__(cls)
        return instance

    def __init__(self, Y, f):
        self.Y = 2*Y #Visto em Benasciutti Tovo sobre one-sided spectral density
        self.f = f

    def moment0(self):
        df = 1

        soma0 = 0

        # Por Dirlik, mi = (1/2pi) integral (w^i)*G(w)dw de -inf até inf

        for i in range(1, len(self.f) - 1):
            if self.f[i] >= 0:
                soma0 += self.Y[i]

        m0 = (self.Y[0] + 2 * soma0 + self.Y[len(self.f) - 1]) * df / 2
        return m0

    def moment1(self):
        soma1 = 0
        df = 1

        for i in range(1, len(self.f) - 1):
            if self.f[i] >= 0:
                soma1 += np.abs(self.Y[i] * self.f[i])

        m1 = (self.Y[0] * self.f[0] + 2 * soma1 + self.Y[len(self.f) - 1] * self.f[len(self.f) - 1]) * df / 2
        return m1

    def moment2(self):
        soma2 = 0
        df = 1

        for i in range(1, len(self.f) - 1):
            if self.f[i] >= 0:
                soma2 += np.abs(self.Y[i] * self.f[i] ** 2)

        m2 = (self.Y[0] * self.f[0] ** 2 + 2 * soma2 + self.Y[len(self.f) - 1] * self.f[len(self.f) - 1] ** 2) * df / 2
        # Sem dividir por 2 o resultado fica mais aceitável
        return m2

    def moment4(self):
        soma4 = 0
        df = 1

        for i in range(1, len(self.f) - 1):
            if self.f[i] >= 0:
                soma4 += np.abs(self.Y[i] * self.f[i] ** 4)

        m4 = (self.Y[0] * self.f[0] ** 4 + 2 * soma4 + self.Y[len(self.f) - 1] * self.f[len(self.f) - 1] ** 4) * df / 2
        return m4

