# Fatiguepy Package

## Available Methods
This package can estimate fatigue life by 6 methods in frequency domain and a method in time domain:

* Rainflow

* Narrow Band
* Rice
* Wirsching-Light
* Tovo-Benasciutti
* <img src="https://render.githubusercontent.com/render/math?math=\alpha_{0.75}">
* Dirlik
* Zhao-Baker

The package requires numpy, math and rainflow

## Instalation

Install directly the package by pip:

```
pip install fatiguepy
```

## Obtaining Power Spectral Density

First, it is necessary to do the calculations of the probability moments.
So, you need Power Spectral Density. To test this package, a random history filtered will be used to get PSD, as seen below

```python
import numpy as np
from scipy import signal

xf = 10
fs = 1024
dt = 1/fs

x = np.arange(0, xf, dt)

signal1 = []
for i in range(len(x)):
    signal1.append(random.random()*300)
signal1 = signal1 - np.mean(signal1)

nyquist = ((len(x)/max(x))/2)

left_pass  = 1.1*100/nyquist
left_stop  = 0.9*100/nyquist
right_pass = 0.9*150/nyquist
right_stop = 1.1*150/nyquist

(N, Wn) = signal.buttord(wp=[left_pass, right_pass],
            ws=[left_stop, right_stop],
            gpass=2, gstop=30, analog=0)

(b, a) = signal.butter(N, Wn, btype='band', analog=0, output='ba')

y = signal.filtfilt(b, a, signal1)

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

E0 = moments.E0()
EP = moments.EP()
gamma = moments.alpha2()
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

The damage intensity (Damage/unit of time) calculated by every method is given by the following equation:

<p align=center>
<img src="https://render.githubusercontent.com/render/math?math=\overline{D} = \nu_pC^{-1}\int_0^\infty s^k p_a(s) ds">
</p>

Where <img src="https://render.githubusercontent.com/render/math?math=p_a"> is PDF of amplitude, k and C is material property, <img src="https://render.githubusercontent.com/render/math?math=\nu_p"> is equivalent to expected number of peaks and s is stress.

## Narrow Band (NB)

For narrow band processes it is reasonable to assume that every peak coincides with a cycle and that, consequently, the amplitudes of the cycles are distributed according to a Rayleigh function.

PDF returns the Probability Density Function of Narrow-Band Method, counting_cycles convert PDF in n [cycles] and loading_spectrum returns the number of cycles having amplitude higher or equal to s.

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
pNB = NB.PDF()

DNB = NB.Damage()
TNB = NB.Life()
TNBh = NB.Lifeh()
```

Damage returns the Damage by NB approach, Life returns the period until failure (in cycles) and Lifeh returns the life in hours.

For the history in study, this method return the respective results:

```
DNB = 7.35e-09 per second
TNB = 1.35e+07 cycles
TNBh = 37772.67 hours
```

## Rice (RC)

Rice model is a mix of Rayleigh and a gaussian distribution.

PDF returns the Probability Density Function of Narrow-Band Method, counting_cycles convert PDF in n [cycles] and loading_spectrum returns the number of cycles having amplitude higher or equal to s.

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

RC = Narrow_Band.RC(k, C, Gyy, f, xf, s)
pRC = RC.PDF()

DRC = RC.Damage()
TRC = RC.Life()
TRCh = RC.Lifeh()
```

For the history in study, this method return the respective results:

```
DNB = 7.09e-09 per second
TNB = 1.41e+07 cycles
TNBh = 39201.35 hours
```

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
For the history in study, this method return the respective results:
```
DWL = 5.76e-09 per second
TWL = 1.74e+07 cycles
TWLh = 48206.58 hours
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
For the history in study, this method return the respective results:
```
DTB = 7.23e-09 per second
TTB = 1.38e+07 cycles
TTBh = 38414.58 hours
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
For the history in study, this method return the respective results:
```
DAL = 7.18e-09 per second
TAL = 1.39e+07 cycles
TALh = 38671.59 hours
```

## Dirlik

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
For the history in study, this method return the respective results:
```
DDK = 7.21e-09 per second
TDK = 1.38e+07 cycles
TDKh = 38539.97 hours
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
For the history in study, this method return the respective results:
```
DZB = 7.15e-09 per second
TZB = 1.39e+07 cycles
TZBh = 38848.97 hours
```

## Rainflow

If you want to calculate rainflow amplitude histogram, you can use the Rainflow module of this package.

### Parameters

<img src="https://render.githubusercontent.com/render/math?math=C"> (*float*):
Constant of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=k"> (*float*):
Slope of SN Curve

<img src="https://render.githubusercontent.com/render/math?math=y"> (*ndarray*):
Stress History y

<img src="https://render.githubusercontent.com/render/math?math=x"> (*ndarray*):
time

The results can be obtained in the same way as the previous methods:

```python
RF = Rainflow.rainflowD(C, k, y, x)

DRF = RF.Damage()
TRF = RF.life()
TRFh = RF.lifeh()
```
For the history in study, this method return the respective results:
```
DRF = 7.05e-09 per second
TRF = 1.42e+07 cycles
TRFh = 39389.11 hours
```

You can compare the result of the methods with counting_cycles() method, present in Narrow_Band, Rice, Dirlik and Zhao_Baker modules as seen below:

```python
import matplotlib.pyplot as plt

S, nRF = RF.rainflow_histogram()

nNB = NB.counting_cycles()
nRC = RC.counting_cycles()
nDK = DK.counting_cycles()
nZB = ZB.counting_cycles()

plt.figure("Comparison in same plot")
plt.bar(S, p,width=round(max(S)*0.015, 1), color='white', edgecolor='black')
plt.plot(s, nNB, linestyle=':',color='blue')
plt.plot(s, nRC, linestyle='-.',color='red')
plt.plot(s, ns, linestyle='-',color='black')
plt.plot(s, nZB, linestyle='--',color='purple')
plt.legend(("Narrow Band", "Rice", "Dirlik", "Zhao-Baker", "Rainflow"))
plt.xlabel(r'Amplitude [MPa]')
plt.ylabel(r'n [cycles]')
plt.grid(True)
plt.show()
```

![Comparison between methods in frequency domain and Rainflow](image/Counting_cycles.png)

It's possible to compare the methods with loading_spectrum as well.

```python
CC, Ss = RF.CumuCycles()

CNB = NB.loading_spectrum()
CRC = NB.loading_spectrum()
CDK = NB.loading_spectrum()
CZB = NB.loading_spectrum()

plt.figure("Comparison between Cumulative Cycles")

plt.semilogx(CC, S, marker='o')
plt.semilogx(CNB, s, marker='D')
plt.semilogx(CRC, s, marker='*')
plt.semilogx(Cs, s, marker='^')
plt.semilogx(CZB, s, marker='s')

plt.legend(("RF", "NB", "RC", "DK", "ZB"))

plt.xlabel('Cumulated Cycles [cycles]')
plt.ylabel(r'S$_a$ [MPa]')

plt.grid(True)
plt.show()
```

![Comparison between methods in frequency domain and Rainflow](image/Cumulative_Cycles.png)

## Relative Error

To compute relative error of any method, the relative_error function, present in all modules of the fatiguepy package, must be used. 

This relative error is in relation to Damage/(unit of second) when type='damage' or in relation to Life when type='cycles'.

Here's an example, calculating error for Zhao-Baker Method:

```python
ZB = Zhao_Baker.ZB(k, C, Gyy, w, xf, s, type='cycles')
psZB = ZB.PDF()

DZB = ZB.Damage()
err = ZB.relative_error(y)
```

When the method parameter is hidden, method="Rainflow" is considered.

If you want to calculate error in relation to the experimental result, do as follows (experimental_value has to be in Damage/(unit of time or unit of life (cycles, s, h, etc))):

```python
ZB = Zhao_Baker.ZB(k, C, Gyy, w, xf, s)
psZB = ZB.PDF()

DZB = ZB.Damage()
Dex = 7.00e-09
err = ZB.relative_error(y, method="Experimental", experimental_value = Dex)
```

Access https://github.com/thiagonds65/fatiguepy to view comparison images
