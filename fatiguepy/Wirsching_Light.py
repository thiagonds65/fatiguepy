import numpy as np

class WL:
    def __new__(cls, *args, **kwargs):
        instance = super(WL, cls).__new__(cls)
        return instance

    def __init__(self, alpha2, k, DNB, xf):
        self.alpha2 = alpha2
        self.k = k
        self.DNB = DNB
        self.xf = xf

    def Damage(self):
        epsylon = np.sqrt(np.abs(1 - self.alpha2 ** 2))

        ak = 0.926 - 0.033 * self.k
        bk = 1.587 * self.k - 2.323

        roWL = ak + (1 - ak) * (1 - epsylon) ** bk

        DWL = roWL * self.DNB

        return DWL

    def Life(self):  # ciclos
        TWL = 1 / WL(self.alpha2, self.k, self.DNB, self.xf).Damage()
        return TWL

    def Lifeh(self):  # horas
        TWLh = self.xf * WL(self.alpha2, self.k, self.DNB, self.xf).Life() / (60 * 60)
        return TWLh
