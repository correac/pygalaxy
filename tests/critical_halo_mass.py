#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as pyplot

z = np.arange(0,6.2,0.2)
Mcritz = []
initial_quess = 12.0
for zi in z:
    Mi = fsolve(Mcrit, initial_quess, args=(zi))
    Mcritz = np.append(Mcritz,Mi)
    initial_guess = Mi

fig = plt.figure()
sub = plt.subplot(1,1,1)
sub.grid(True)
plt.plot(z,Mcritz)
plt.xlabel('$z$')
plt.ylabel('$\log_{10}M_{200}[M_{\odot}]$')
#plt.axis([1e4,1e9,1e-24,1e-21])
