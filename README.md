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
>>>import numpy as np

>>>x = np.linspace(0, 1000, 100000)
>>>y = 100 * np.sin(2 * np.pi * 5 * x) + 70 * np.cos(w2 = 2 * np.pi * 10 * x) + 90 * np.sin(2 * np.pi * 20 * x)
>>>Y = np.fft.fft(y)
>>>Y = 2.0 * np.abs(Y / len(Y))
>>>f = np.fft.fftfreq(len(Y)) * n / xf * 2 * np.pi    
>>>mask = f > 0
```

So, just use the module present in the fatiguepy package. Function moment0 to moment4 returns respective probability moment, E0 returns the expected positive zero-crossing rate, EP returns the expected peak occurrency frequency and alpha2 returns spectral width parameter.

```
>>>from fatiguepy import *
>>>moments = prob_moment.Probability_Moment(Y[mask], f[mask])
>>>
>>>m0 = moments.moment0()
>>>m1 = moments.moment1()
>>>m2 = moments.moment2()
>>>m4 = moments.moment4()
>>>m75 = moments.moment0dot75()
>>>m15 = moments.moment1dot5()
>>>
>>>E0num = moments.E0()
>>>EPnum = moments.EP()
>>>gammanum = moments.alpha2()
```

## Narrow Band (NB)

For narrow band processes it is reasonable to assume that every peak coincides with a cycle and that, consequently, the amplitudes of the cycles are distributed according to a Rayleigh function.

PDPeaks returns the Probability Density Function of Narrow-Band Method.

```
>>>si = 1
>>>sf = 200
>>>ds = 1
>>>s = np.arange(si, sf, ds)
>>>
>>>NB = Narrow_Band.NB(k, C, Y[mask], f[mask], xf)
>>>pp = NB.PDPeaks(s)
```

The damage calculated by Narrow Band method is given by Palmgren-Miner applied to Probability Density Function or given by the following equation:

<p align=center>
<img src="https://render.githubusercontent.com/render/math?math=\overline{D}_{NB} = \nu_0C^{-1}\left(\sqrt{2m_0}\right)^k\Gamma\left(1 %2B \frac{k}{2}\right)">
</p>
Damage returns the Damage by NB approach, Life returns the period (in cycles) and Lifh returns the life in hours.

```
>>>NB.Damage()
0.002705610306446634
>>>NB.Life()
369.6023768157997
>>>NB.Lifeh()
102.66732689327769
>>>
>>>nt = xf * EPnum
>>>
>>>npi = pp * ds * nt
>>>
>>>Nfi = 0.5 * (s/sigmaf) ** (1/b)
>>>
>>>Dp = 0
>>>
>>>for i in range(len(s)):
...    Dp += npi[i] / Npfi[i]
>>>Dp
0.0027315862098964768
>>>(xf/Dp)/3600
101.69101629353487
```

## Wirsching-Light (WL)

To this method, Wirsching and Light considered an width parameter to correct Narrow-Band approximation with an empirical factor. It can be done with the fatiguepy package as follows:

```
>>>WL = Wirsching_Light.WL(k, C, Y[mask], f[mask], xf)
>>>
>>>WL.Damage()
0.0018584040811990177
>>>WL.Life()
538.096106286429
>>>WL.Lifeh()
149.47114063511916
```

## Dirlik (DK)

This method has long been considered to be one of the best and has already been subject to modifications, e.g., for the inclusion of the temperature effect.

The functions for this method are analogous to the NB functions:

```
>>>DK = Dirlik.DK(k, C, Y[mask], f[mask], xf)
>>>
>>>ps = DK.PDF(s)

>>>DK.Damage()
0.04298062994937626
>>>DK.Life()
23.266294634997834
>>>DK.Lifeh()
6.462859620832732
>>>nt = xf * EPnum
>>>
>>>ni = ps * ds * nt
>>>
>>>Nfi = 10 ** ((2.967079 - np.log10(s)) / 0.138)
>>>
>>>D = 0
>>>
>>>for i in range(len(s)):
...    D += ni[i] / Nfi[i]
>>>D
0.10769442833870352
>>>(xf / D) / 3600
2.579314288239267
```
