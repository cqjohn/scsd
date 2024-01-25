#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 16:46:41 2024

@author: genah
"""
#import libraries 
import numpy as np 
import matplotlib.pyplot as plt 
import random 
from scipy.fft import fft, fftfreq, ifft


#generate x axis values in both domains 
dt = 0.001
t = np.arange(0,100,dt)
n = len(t)
freq = (1/(dt*n)) * np.arange(0.00000000001,n)

#generate random 1/f noise
def pink_noise(f,alph=1):
    
    X = [random.uniform(0, 2*np.pi) for _ in range(0,len(f))]
    Y = [random.uniform(0, 2*np.pi) for _ in range(0,len(f))]
    mag = [np.sqrt(X[i]**2 + Y[i]**2) for i in range(0,len(f))]
    u =  mag * (f**(-alph)) 
    return u 

#calcualte PSD 
u = pink_noise(freq)


plt.figure(1)
plt.plot(freq,u)
plt.xlim(10**(-2))
plt.ylim(10**(-8),10**6)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Fequenecy (Hz)")
plt.ylabel("Spectral Density $(\mu V^2/Hz)$")
plt.show()

#%%
cgam = ifft(u) #inverse fourier transform 
#plotc = (cgam * np.conjugate(cgam))/n  
plt.figure(2)
plt.plot(t,cgam)
plt.xlabel("Time (s)")
plt.ylabel("Current")
plt.show()

#%%
#can i get it back 
original_signal= fft(cgam)
plt.figure(3)
plt.plot(freq, original_signal)
plt.xlim(10**(-2))
plt.ylim(10**(-8),10**6)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Fequenecy (Hz)")
plt.ylabel("Spectral Density $(\mu V^2/Hz)$")
plt.show()


#questions psd has a sqaure value do I need to squre 