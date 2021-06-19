from collections import defaultdict
import fatpack
import numpy as np

class rainflowD:

    def __init__(self, C, k, y, x, nbins=50):
        self.C = C
        self.k = k
        self.b = -1/k
        self.A = C**(1/k)
        self.sigmaf = self.A/(2**self.b)
        self.y = y
        self.x = x
        self.nbins = nbins

    def rainflow_histogram(self):
        # Preciso otimizar essa função inteira pra conseguir gerar o sm e o r do histórico. Tentar aprender a usar numba depois
        r, sm = fatpack.find_rainflow_ranges(self.y, return_means=True, k=self.nbins)

        # Next create the bins to divide the stress ranges and means into

        bins_r = np.linspace(min(r), max(r), self.nbins)
        bins_sm = np.linspace(min(sm), max(sm), self.nbins)

        # and establish the data array. Note that other datavectors vectors, e.g. 
        # Smin and Smax, may also be used to create other data arrays and
        # resulting rainflow matrices.

        data_array = np.array([sm, r]).T

        # Finally, establish the rainflow matrix from the data array and the
        # specified row and column bins.

        n = fatpack.find_rainflow_matrix(data_array, bins_sm, bins_r)
        n = n.sum(axis=0)
        r = bins_r[0:-1]
        S = (r/2)/(1-np.mean(self.y)/self.sigmaf)
        # r = []
        # n = []

        # for i in range(len(cc)):
        #     r.append(cc[i][0]/2)
        #     n.append(cc[i][1])
        #     # print(r[i], n[i])

        # O código abaixo acho que está certo, porém preciso dar um jeito de fazê-lo ser mais rápido para obter os ranges e means de rainflow

        rangemax = max(S) - min(S)

        nclass = self.nbins
        tns = sum(n)

        bw = (min(S) + (rangemax/nclass)*(nclass+1)) - (min(S) + (rangemax/nclass)*(nclass))

        p = []
        appendp = p.append

        summ = 0
        for i in range(len(n)):
            appendp((n[i]/tns)/bw)
            summ += p[i]*bw
        
        #Finish here!!!!

        return S, n, p

    def CumuCycles(self):
        S, n, p = self.rainflow_histogram()
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

        S, n, p = self.rainflow_histogram()
        DRF = 0
        for i in range(1,len(n)):
            Nf = 0.5*(S[i]/self.sigmaf)**(1/self.b)
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

