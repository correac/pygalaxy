#!/usr/bin/env ipython
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pyplot
from wiersma_cooling_table import compute_net_cooling_normal_opt

x = np.arange(4,9.1,0.02)
x = 10**x
y_0 = []
y_2 = []
y_4 = []
y_6 = []
for T in x:
    y_0 = np.append(y_0,compute_net_cooling_normal_opt(redshift=3.017,density=10**(0),temperature=T,metallicity=1.0))
    y_2 = np.append(y_2,compute_net_cooling_normal_opt(redshift=3.017,density=10**(-2),temperature=T,metallicity=1.0))
    y_4 = np.append(y_4,compute_net_cooling_normal_opt(redshift=3.017,density=10**(-4),temperature=T,metallicity=1.0))
    y_6 = np.append(y_6,compute_net_cooling_normal_opt(redshift=3.017,density=10**(-6),temperature=T,metallicity=1.0))

fig = plt.figure()
sub = plt.subplot(1,1,1)
sub.grid(True)
plt.plot(x,np.abs(y_0),label='$n_{H}=10^{0}$ cm$^{-3}$')
plt.plot(x,np.abs(y_2),label='$n_{H}=10^{-2}$ cm$^{-3}$')
plt.plot(x,np.abs(y_4),label='$n_{H}=10^{-4}$ cm$^{-3}$')
plt.plot(x,np.abs(y_6),label='$n_{H}=10^{-6}$ cm$^{-3}$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('$|\Lambda/n_{H}^{2}|$ [erg cm$^{3}$/s]')
plt.xlabel('T [K]')
plt.legend()
plt.axis([1e4,1e9,1e-24,1e-21])
