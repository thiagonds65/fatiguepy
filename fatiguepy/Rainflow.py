import rainflow
import matplotlib.pyplot as plt
import numpy as np
from . import prob_moment

class rainflowD:
    def __new__(cls, *args, **kwargs):
        instance = super(rainflowD, cls).__new__(cls)
        return instance

    def __init__(self, C, k, y, xf):
        self.C = C
        self.k = k
        self.b = -1/k
        self.A = C**(1/k)
        self.y = y
        self.xf = xf

    def PDF(self, Y, f):
        cc = rainflow.count_cycles(self.y)

        '''
        nn = np.zeros(len(cc))
        rng = np.zeros(len(cc))
        nnn = np.zeros(len(cc))
        rngrou = np.zeros(len(cc))
    
        for j in range(len(cc)):
            nn[j] = cc[j][1]
            rng[j] = cc[j][0]
            for i in range(1, len(cc)):
                if round(rng[i - 1], 4) == round(rng[j - 1], 4):
                    rngrou[j] = round(rng[i], 4)
                    nnn[j] += nn[i]
        '''

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
        
        
        for i in range(len(s) - 1):
            for j in range(len(s) - 1):
                if s[j] > s[j+1]:
                    auxs = s[j]
                    s[j] = s[j+1]
                    s[j+1] = auxs

                    auxn = n[j]
                    n[j] = n[j+1]
                    n[j+1] = auxn
        
        '''
        p = []
        ss = []
        c = 1
        for i in range(0, len(s)-1):
            if round(s[i], 4) == round(s[i+1], 4):
                c +=1
                
                sss = s[i]
            else:
                p.append(c/len(s))
                ss.append(sss)
                c = 1
        '''
        
        
        p = np.zeros(len(s))
        ds = s[1] - s[0]
        
        EP = prob_moment.Probability_Moment(Y, f).EP()
            
        for i in range(len(s)):
            p[i] = n[i]/(ds*EP*self.xf)
        
        ds = s[1] - s[0]
        DRF = 0
        for i in range(1,len(p)):
            DRF += (EP*(self.C**(-1))*(s[i]**self.k)*p[i]*ds)
        
        
        print(f"Dano Rainflow PDF: {DRF}")

        '''
        plt.figure(1)
        plt.plot(s, p)

        plt.xlabel(r'\Delta{S} [MPa]')
        plt.ylabel(r'PDF [MPa$^{-1}$]')

        plt.grid(True)
        plt.show()
        '''
        
        return s, p

    def Damage(self):
        cc = rainflow.count_cycles(self.y)

        vl = []
        pk = []
        cnt = []
        sa = []
        sm = []
        j = 0
        
            
        for valley, peak, count in rainflow.extract_cycles(self.y): 
            #print(valley, peak, count)
            vl.append(valley)
            pk.append(peak)
            cnt.append(count)
            sa.append((pk[j]-vl[j])/2) 
            sm.append((pk[j]+vl[j])/2)
                
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


        DRF = DRF/self.xf
        
        return DRF
    
    def Lifes(self):
        TRFs =  1/self.Damage()
        return TRFs
    
    def Lifeh(self):
        TRFh = self.Lifes()/3600
        return TRFh
    
    def Life(self):
        TRF = self.Lifes()/self.xf
        return TRF
