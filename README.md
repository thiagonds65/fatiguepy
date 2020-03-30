# Fatiguepy Package

## Available Methods
This package can estimate fatigue life by 4 methods:

* Narrow Band
* Wirsching-Light
* Dirlik
* Zhao-Baker

## Probability Moments

First, it is necessary to do the calculations of the probability moments.
So, you need Power Spectral Density. To test this package, sum of sinusoid will be used to get PSD, as seen below

```
import numpy as np
from fatiguepy import *

w1 = 2 * np.pi * 5
w2 = 2 * np.pi * 10
w3 = 2 * np.pi * 20
xf = 1000
fa = 1/xf
n = 100 * xf

x = np.linspace(0, xf, n)
y = 100 * np.sin(w1 * x) + 70 * np.cos(w2 * x) + 90 * np.sin(w3 * x)

Y = np.fft.fft(y)

Y = 2.0 * np.abs(Y / len(Y))

f = np.fft.fftfreq(len(Y)) * n / xf * 2 * np.pi    

mask = f > 0
```
So, just use the module present in the fatiguepy package

```
moments = prob_moment.Probability_Moment(Y[mask], f[mask])

m0 = moments.moment0()
m1 = moments.moment1()
m2 = moments.moment2()
m4 = moments.moment4()
m75 = moments.moment0dot75()
m15 = moments.moment1dot5()

E0num = moments.E0()

EPnum = moments.EP()

gammanum = moments.alpha2()
```

## Narrow Band (NB)

For narrow band processes it is reasonable to assume that every peak coincides with a cycle and that, consequently, the amplitudes of the cycles are distributed according to a Rayleigh function.

PDF is obtained through the following lines of code:

```
si = 1
sf = 200
ds = 1
s = np.arange(si, sf, ds)

NB = Narrow_Band.NB(k, C, Y[mask], f[mask], xf)
pp = NB.PDPeaks(s)
```

The damage calculated by Narrow Band method is given by Palmgren-Miner applied to graph or given by the following equation:

<p align=center>
<img src="https://render.githubusercontent.com/render/math?math=\overline{D}_{NB} = \nu_0C^{-1}\left(\sqrt{2m_0}\right)^k\Gamma\left(1 - \frac{k}{2}\right)">
</p>
Following code shows how this is done.

```
DNB = NB.Damage()
TNB = NB.Life()
TNBh = NB.Lifeh()

nt = xf * EPnum

npi = pp * ds * nt

Nfi = 0.5 * (s/sigmaf) ** (1/b)

Dp = 0

for i in range(len(s)):
    Dp += npi[i] / Npfi[i]
```