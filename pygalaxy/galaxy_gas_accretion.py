#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import numpy as np
import scipy.integrate as integrate

def density_model(r,M200,z):
    function = lambda x, a, b, c: a + b*x + c*x**2
    if z<=2:
        gamma_1 = function(z,-0.842,-1.588, 0.297)
        gamma_2 = function(z,-0.629,0.523,-0.114)
        gamma_3 = function(z,-0.282,-0.481,0.361)
        beta_1 = function(z,0.323,-0.819,0.148)
        beta_2 = function(z,-1.035,0.472,-0.119)
        beta_3 = function(z,0.052,-0.724,0.388)
        gamma = function(M200-12.0,gamma_1,gamma_2,gamma_3)
        beta = function(M200-12.0,beta_1,beta_2,beta_3)
    else:
        gamma = function(M200-12.0,-2.83,-0.039,0.2)
        beta = function(M200-12.0,-0.723,-0.567,0.156)
    
    rho_function = r**gamma * 10**(beta * (np.log10(r))**2) #r in units of R200
    return rho_function

def integral(x,M200,z):
    f = 4. * np.pi * density_model(x,M200,z) * x**2
    return f

def density_profile_bestfit_to_EAGLE(r,M200,z):
    
    Msun = 1.98847e33 #Msun->gr
    f_baryon = 0.04825/0.307
    R200 = 10**M200 * Msun / (4. * np.pi * 200. * rhocrit(z) / 3.) #=R200^3 [cm^3]
    #rho_norm = 10**M200 * Msun * f_baryon / R200
    rho_norm = 10**Mhot(M200,z) * Msun / R200
    norm = integrate.quad(integral, 0., 1., args=(M200,z))[0]
    rho_norm /= norm
    
    
    rho_hot_gas = rho_norm * density_model(r,M200,z) #r in units of R200, rho [gr/cm^3]
    return rho_hot_gas

def density_profile_isothermal(r,M200,z):
    Msun = 1.98847e33 #Msun->gr
    f_baryon = 0.04825/0.307
    R200 = 10**M200 * Msun / (4. * np.pi * 200. * rhocrit(z) / 3.) #=R200^3 [cm^3]
    rho_norm = 0.9 * f_baryon * 10**M200 * Msun / (4. * np.pi * R200)
    rho_hot_gas = rho_norm / r**2  #r in units of R200, rho [gr/cm^3]
    return rho_hot_gas


def Gamma_cool_with_radial_dependence(M200,z,r,Thot=1e6,Zhot=0.1):
    """
        Function that calculates the radiative cooling rate as defined in
        eq. (22) of Correa et al. (2018), MNRAS, 473, 1, 538.
        
        Parameters
        ----------
        M200 : float
        Logarithmic halo mass at redshift zero. Note, the halo is in
        units of M200 (mass enclosed within the virial radius, R200,
        defined as the radius at which the mean density is 200 times
        the critical).
        
        z : float / array
        Redshift to solve acc_rate.
        
        f_baryon : float, optional.
        Universal baryon fraction.
        
        Thot : float, optional.
        Gas temperature. Default halo's virial temperature.
        
        rho_hot : float, optional.
        Gas density. Default 10^0.6 times rho_crit(z).
        
        Zhot : float, optional.
        Gas metallicity. Default 0.1 times solar metallicity.
        
        Returns
        -------
        Gamma_cool : float / array
        Radiative cooling rate [erg/s] at redshift 'z'.
        
        Notes
        -----
        This functions reads the cooling table from Wiersma et al. (2009).
        """
    Msun = 1.98847e33 #Msun->gr
    Thot = Tvir(M200,z) #hot gas temperature [K]
    rho_hot = density_profile_bestfit_to_EAGLE(r,M200,z) #density of hot gas (gr/cm^3)
    
    Gamma = 10**Mhot(M200,z) * Msun #units gr
    Gamma *= Lambda(z, rho_hot[0], Thot, Zhot) #units gr erg cm^-3 s^-1
    Gamma /= rho_hot #units erg/s
    return Gamma

def equaling_heating_cooling_rates(x,M200,z):
    f = Gamma_heat(M200,z)-Gamma_cool_with_radial_dependence(M200,z,x)
    return f

def cooling_radius(M200,z):
    initial_quess = 0.5
    rcool = fsolve(equaling_heating_cooling_rates, initial_quess, args=(M200,z))
    return rcool

def galaxy_gas_accretion(M200,z):
    
    rcool = cooling_radius(M200,z) # in units of R200
    epsilon = 0.3
    f_halo_acc_hot = fhotacc(M200,z)
    f_halo_acc_cold = 1.0-f_halo_acc_hot
    dMdt_gas_halo = M200accr(M200,z,commah_flag=False) #Msun/yr
    dMdt_gas_galaxy = epsilon * (f_halo_acc_hot * rcool+f_halo_acc_cold) * dMdt_gas_halo #Msun/yr
    return dMdt_gas_galaxy

