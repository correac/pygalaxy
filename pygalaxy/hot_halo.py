#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import numpy as np
import commah
from wiersma_cooling_table import Lambda

def fhot(M200,z):
    x = M200-12.
    z_tilde = np.log10(1.+z)
    if z<=2:
        alpha = -0.79+0.31*z_tilde-0.96*z_tilde**2
        beta = 0.52-0.57*z_tilde+0.85*z_tilde**2
        gamma = -0.05
    else:
        alpha = -0.38-1.56*z_tilde+1.17*z_tilde**2
        beta = 0.12+0.94*z_tilde-0.55*z_tilde**2
        gamma = -0.05
    return alpha+beta*x+gamma*x**2 #log10_fhot

def Mhot(M200,z,f_baryon=None):
    log10_fhot = fhot(M200,z)
    if f_baryon==None:f_baryon = 0.04825/0.307
    Mhot = log10_fhot+np.log10(f_baryon)+M200
    return Mhot

def fhotacc(M200,z):
    z_tilde = np.log10(1.+z)
    if z<2:
        M_half=-0.15+0.22*z+0.07*z**2
        alpha = -1.86*10**(-1.26*z_tilde+1.29*z_tilde**2)
    if z>=2 & z<4:
        M_half=-0.25+0.53*z-0.07*z**2
        alpha = -0.46*10**(0.81*z_tilde-0.42*z_tilde**2)
    if z>=4:
        M_half = 0.72+0.01*z
        alpha = -1.07
    M_half += 12.0
    x = M200-M_half
    x = 10**x
    fhotacc = 1./(1.+x**alpha)
    return fhotacc

def rhocrit(z,Omega_m=None,h=None):
    if Omega_m==None:Omega_m=0.307
    Omega_l = 1.0-Omega_m
    if h==None:h=0.6777
    H0 = h * 100. #km/s/Mpc
    G = 4.30091e-3 #pc Msun^-1 (km/s)**2
    G *= 1e-3 #pc->Mpc (Mpc Msun^-1 (km/s)**2
    rho0 = 3*H0**2/(8.*np.pi*G) #Msun Mpc^-3
    return rho0 * (Omega_m*(1+z)**3+Omega_l)

def M200accr(M200,z,**cosmo):
    output = commah.run(cosmo, zi=0, M200, z)
    return output['dMdt']

def Tvir(M200,z):
    return 10**5.3 * (10**(M200-12.0))**(2./3.) * (1.+z)

def Gamma_heat(M200,z,**cosmo):
    Msun = 1.98847e33 #Msun->gr
    yr_to_sec = 3.154e7
    mu = 0.59
    k_b = 1.3806485e-23 #m^2 kg s^-2 K^-1
    mp = 1.6726219e-27 #kg
    cte = k_b / (mu * mp) #(m/s)^2 K^-1
    cte *= 1e2 #m->cm, (cm/s)^2 K^-1
    Gamma = (3./2.) * cte * Tvir(M200,z) * f_baryon #(cm/s)^2
    Gamma *= M200accr(M200,z,**cosmo) #(cm/s)^2 Msun/yr
    Gamma *= Msun/yr_to_sec #erg/s
    Gamma *= ((2./3.)*fhot(M200,z) +fhotacc(M200,z))
    return Gamma

def Lambda(z, rho, T, metallicity, Xh=0.752, Xhe = 0.248):
    net_cool = compute_net_cooling_normal_opt(z, rho, T, metallicity) #return Lambda_net/n_H^2 [erg cm^3 s^-1]#
    mu = 0.59
    mp = 1.7e-24     # proton mass: gr
    n_gas = rho/(mu * mp) # log gas number density [cm^-3]
    n_H = n_gas*Xh
    net_cool *= n_H^2 #[erg cm^-3 s^-1]
    return net_cool

def Gamma_cool(M200,z,f_baryon=None,Thot=None,rho_hot=None,Zhot=None):
    if f_baryon==None:f_baryon=0.04825/0.307
    if Thot==None:Thot=Tvir(M200,z)
    if rho_hot==None:rho_hot=10**0.6*rhocrit(z)
    if Zhot==None:Zhot=0.1 #Solar metallicity
    Gamma = fhot(M200,z) * f_baryon * M200 * Lambda(z, rho, Thot, Zhot) / rho_hot
    return Gamma

