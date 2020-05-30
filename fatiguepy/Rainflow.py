import rainflow
import matplotlib.pyplot as plt

class rainflowD:
    def __new__(cls, *args, **kwargs):
        instance = super(rainflowD, cls).__new__(cls)
        return instance

    def __init__(self, A, b, y, xf):
        self.A = A
        self.b = b
        self.y = y
        self.xf = xf
    def DRF(self):
        cc = rainflow.count_cycles(self.y)

        vl = []
        pk = []
        cnt = []
        sa = []
        sm = []
        j = 0
        with open('rainflow.txt', 'w') as file:
            file.write("      Valley                Peak         Count           Sa                 Sm\n")
            for valley, peak, count in rainflow.extract_cycles(self.y): 
                #print(valley, peak, count)
                vl.append(valley)
                pk.append(peak)
                cnt.append(count)
                sa.append((pk[j]-vl[j])/2) 
                sm.append((pk[j]+vl[j])/2)
                
                file.write(str(vl[j])+"    ")
                file.write(str(pk[j])+"   ")
                file.write(str(cnt[j])+"    ")
                file.write(str(sa[j])+"    ")
                file.write(str(sm[j])+" ")
                file.write("\n")
                
                j+=1

        DRF = 0
        Nf = []
        s = []
        n = []
        sigmaf = self.A/(2**self.b)
        for i in range(len(sa)):
            s.append(sa[i]/(1-(sm[i]/sigmaf)))
            n.append(cnt[i])
            Nf.append((s[i]/self.A)**(1/self.b))
            DRF += n[i] / Nf[i]
            # print(f"{n[i]} e {Nf[i]}")

        DRF = DRF/self.xf
        return DRF