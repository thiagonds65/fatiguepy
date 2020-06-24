# Fatiguepy Package

## Available Methods
This package can estimate fatigue life by 5 methods:

* Narrow Band
* Wirsching-Light
* Tovo-Benasciutti
* <img src="https://render.githubusercontent.com/render/math?math=\alpha_{0.75}">
* Dirlik (Rainflow Range and Ordinary Range)
* Zhao-Baker

The package requires numpy, math and rainflow

## Instalation

Install directly the package by pip:

```
pip install fatiguepy
```

## Obtaining Power Spectral Density

First, it is necessary to do the calculations of the probability moments.
So, you need Power Spectral Density. To test this package, sum of sinusoid will be used to get PSD, as seen below

```python
import numpy as np
from scipy import signal

xf = 10
fs = 1024
dt = 1/fs

x = np.arange(0, xf, dt)

y = np.zeros(len(x))

for i in range(100):
    y += random.randint(0,10) * np.sin(2*np.pi*random.randint(0, 100) * x)

window = signal.hann(len(y), False)

f, Gyy = signal.welch(y, fs, return_onesided=True, window=window, average='median')
```

## Probability Moments

Once the PSD and frequency are obtained, just use the module present in the fatiguepy package. Function moment0 to moment4 returns respective probability moment, E0 returns the expected positive zero-crossing rate, EP returns the expected peak occurrency frequency and alpha2 returns spectral width parameter.


### Parameters

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
array of sample frequencies


```python
from fatiguepy import *
moments = prob_moment.Probability_Moment(Gyy, f)

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

Steel SAE 1015 was considered, so Python can perform the calculations.

```python
b = -0.138
sigmaf = 1020
A = (2**b)*sigmaf
k = -1/b
C = A ** k
```

## Damage

The damage (Damage/unit of time) calculated by every method is given by Palmgren-Miner Rule applied to Probability Density Function or given by the following equation:

<p align=center>
<img src="https://render.githubusercontent.com/render/math?math=\overline{D} = \nu_pC^{-1}\int_0^\infty s^k p_a(s) ds">
</p>

Where <img src="https://render.githubusercontent.com/render/math?math=p_a"> is PDF of amplitude, k and C is material property, <img src="https://render.githubusercontent.com/render/math?math=\nu_p"> is equivalent to expected number of peaks and s is stress.

## Narrow Band (NB)

For narrow band processes it is reasonable to assume that every peak coincides with a cycle and that, consequently, the amplitudes of the cycles are distributed according to a Rayleigh function.

PDPeaks returns the Probability Density Function of Narrow-Band Method.

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

NB = Narrow_Band.NB(k, C, Gyy, f, xf, s)
pp = NB.PDPeaks()

DNB = NB.Damage()
TNB = NB.Life()
TNBh = NB.Lifeh()
```

Damage returns the Damage by NB approach, Life returns the period (in cycles) and Lifeh returns the life in hours.

## Wirsching-Light (WL)

To this method, Wirsching and Light considered an width parameter to correct Narrow-Band approximation with an empirical factor. It can be done with the fatiguepy package as follows:

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

WL = Wirsching_Light.WL(k, C, Gyy, f, xf, s)
DWL = WL.Damage()
TWL = WL.Life()
TWLh = WL.Lifeh()
```

## Tovo-Benasciutti (TB)

To this method, Tovo and Benasciutti proposed an approach where the
fatigue life is calculated as a linear combination of the upper and lower fatigue-
damage intensity limits. It can be done with the fatiguepy package as follows:

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

TB = Tovo_Benasciutti.TB(k, C, Gyy, f, xf, s)
DTB = TB.Damage()
TTB = TB.Life()
TTBh = TB.Lifeh()
```

## <img src="https://render.githubusercontent.com/render/math?math=\alpha_{0.75}"> method (AL)

This method is a correction method based in a spectral parameter <img src="https://render.githubusercontent.com/render/math?math=\alpha_{0.75}">, and it's can be done as follows:

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=xf"> (*float*):
Observation period

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

AL = alpha075.AL(k, C, Gyy, f, xf, s)
DAL = AL.Damage()
TAL = AL.Life()
TALh = AL.Lifeh()
```

## Dirlik *Ordinary Range Half Cycle* (OR)

The ordinary range behaves in small ranges like an exponential decrease close to origin. The later part of the densities features a Rayleigh Function.

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


This method works as seen below:

```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

DK = Dirlik.DK(k, C, Gyy, f, xf, s)

psOR = DK.PDFOR()

DOR = DK.DamageOR()
TOR = DK.LifeOR()
TORh = DK.LifehOR()
```

## Dirlik *Rainflow Range Half Cycle* (RR)

This method has long been considered to be one of the best and has already been subject to modifications, e.g., for the inclusion of the temperature effect.

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


The functions for this method are analogous to the NB functions:

```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

DK = Dirlik.DK(k, C, Gyy, f, xf, s)

ps = DK.PDF()

DDK = DK.Damage()
TDK = DK.Life()
TDKh = DK.Lifeh()
```

## Zhao-Baker (ZB)

This method combined theoretical assumptions and simulation results to give the linear combination of Weibull and Rayleigh Probability Density Function.

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=G_{yy}"> (*ndarray*):
Power Spectral Density or Power Spectrum of Stress History y

<img src="https://render.githubusercontent.com/render/math?math=f"> (*ndarray*):
Array of sample frequencies

<img src="https://render.githubusercontent.com/render/math?math=s"> (*ndarray*):
Array of sample stresses


The results can be obtained in the same way as the previous methods:

```python
si = 0.0
sf = abs(max(y)-min(y))
ds = sf/128
s = np.arange(si, sf, ds)

ZB = Zhao_Baker.ZB(k, C, Gyy, w, xf, s)
psZB = ZB.PDF()

DZB = ZB.Damage()
TZB = ZB.Life()
TZBh = ZB.Lifeh()
```

## Relative Error

To compute relative error of any method, the relative_error function, present in all modules of the fatiguepy package, must be used, with the exception of the Dirlik module, which has a difference between the Rainflow Range and Ordinary Range methods. In these cases, you must use the relative_errorRR() method, for Rainflow Range, and relative_errorOR(), for Ordinary Range. 

This relative error is in relation to Damage/(unit of second).

Here's an example, calculating error for Zhao-Baker Method:

```python
ZB = Zhao_Baker.ZB(k, C, Gyy, w, xf, s)
psZB = ZB.PDF()

DZB = ZB.Damage()
err = ZB.relative_error(y)
```

When the method parameter is hidden, method="Rainflow" is considered.

If you want to calculate error in relation to the experimental result, do as follows (Dexperimental has to be in Damage/(unit of time)):
```python
ZB = Zhao_Baker.ZB(k, C, Gyy, w, xf, s)
psZB = ZB.PDF()

DZB = ZB.Damage()
Dex = 0.61
err = ZB.relative_error(y, method="Experimental", Dexperimental = Dex)
```
