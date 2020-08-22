import rainflow
import matplotlib.pyplot as plt
import numpy as np
from . import prob_moment

class rainflowD:


    def __init__(self, C, k, y, x):
        self.C = C
        self.k = k
        self.b = -1/k
        self.A = C**(1/k)
        self.y = y
        self.x = x

    def rainflow_histogram(self):
        cc = rainflow.count_cycles(self.y, nbins=50)

        r = []
        n = []
        for i in range(len(cc)):
            r.append(cc[i][0]/2)
            n.append(cc[i][1])

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
        '''

        # Begins here!!!!
        tns = sum(n)

        rangemax = max(r) - min(r)

        nclass = 50

        rangesm = []
        auxs = []
        smean = []

        nmatrix = []
        S = []
        sigmaf = self.A/(2**self.b)

        for j in range(nclass):
            sumn = 0
            for i in range(len(r)):
                if(r[i] >= min(r) + (rangemax/nclass)*(j)  and r[i] < min(r) + (rangemax/nclass)*(j+1)):
                    auxs.append(r[i])
                    sumn += n[i]

                if (r[i] == max(r) and j == nclass-1):
                    auxs.append(r[i])
                    sumn +=n[i]
            
            nmatrix.append(sumn)

            smean.append((min(r) + (rangemax/nclass)*(j) + min(r) + (rangemax/nclass)*(j+1))/2)
            rangesm.append(auxs[:])
            auxs.clear()

            S.append(smean[j]/(1-(self.y.mean()/sigmaf)))



        bw = (min(r) + (rangemax/nclass)*(j+1)) - (min(r) + (rangemax/nclass)*(j))

        p = []
        summ = 0

        for i in range(len(nmatrix)):
            p.append((nmatrix[i]/tns)/bw)
        
        #Finish here!!!!

        return S, nmatrix

    def CumuCycles(self):
        S, n = self.rainflow_histogram()
        CumuC = np.zeros(len(S))

        for i in range(len(S)):
            for j in range(i, len(S)):
                CumuC[i] += n[j]
        
        return CumuC, S


    def Damage(self):
        '''
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


            DRF = DRF/max(self.x)
        '''

        S, n = self.rainflow_histogram()
        DRF = 0
        for i in range(1,len(n)):
            Nf = self.C/(S[i]**self.k)
            DRF += n[i]/Nf 
        
        return DRF/max(self.x)

    def Lifes(self):
        TRFs =  1/self.Damage()
        return TRFs
    
    def Lifeh(self):
        TRFh = self.Lifes()/3600
        return TRFh
    
    def Life(self):
        TRF = self.Lifes()/max(self.x)
        return TRF
