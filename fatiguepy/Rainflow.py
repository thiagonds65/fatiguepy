import rainflow

class rainflowD:
    def __new__(cls, *args, **kwargs):
        instance = super(rainflowD, cls).__new__(cls)
        return instance

    def __init__(self, A, b, y):
        self.A = A
        self.b = b
        self.y = y
        
    def DRF(self):
        cc = rainflow.count_cycles(self.y)

        DRF = 0
        Nf = []
        s = []
        n = []
        for i in range(len(cc)):
            s.append(rainflow.count_cycles(self.y)[i][0])
            n.append(rainflow.count_cycles(self.y)[i][1])
            Nf.append((0.5*s[i]/self.A)**(1/self.b))
            DRF += n[i] / Nf[i]
            # print(f"{n[i]} e {Nf[i]}")

        return DRF