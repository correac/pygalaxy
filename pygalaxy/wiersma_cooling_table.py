#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

import numpy as np
import os
import hdf5
from scipy import interpolate

def read_option(option):
    dir = './wiersma09_coolingtables/'+option+'_option_hdf5_files/'
    for entry in os.listdir(dir):
        if os.path.isfile(os.path.join(dir, entry)):
            file = entry
    return dir+file

def look_for_table_name(redshift):
    dir = './wiersma09_coolingtables/normal_option_hdf5_files/'
    z = []
    file = []
    for entry in os.listdir(dir):
        if os.path.isfile(os.path.join(dir, entry)):
            file = np.append(file,entry)
            a,b = entry.split("_")
            a,b = b.split(".hdf5")
            z = np.append(z,float(a))
    idx = np.argsort(z)
    z = z[idx]
    file = file[idx]
    dz = z-redshift
    min_z = np.min(np.abs(dz))
    select_file = np.where(np.abs(dz)==min_z)[0]
    output = dir+file[select_file[0]]
    if dz[select_file]<0:
        output2 = dir+file[select_file[0]+1]
        dz = np.array([z[select_file[0]],z[select_file[0]+1]])
    else:
        output2 = dir+file[select_file[0]-1]
        dz = np.array([z[select_file[0]],z[select_file[0]-1]])
    return output, output2, dz

def read_wiersma_cooling_table(filename):
    with h5py.File(filename, "r") as hf:
        data = {}
        data['H_bins'] = hf["/Total_Metals/Hydrogen_density_bins"][:]
        data['Temp_bins'] = hf["/Total_Metals/Temperature_bins"][:]
        data['total_metals_net_cooling'] = hf["/Total_Metals/Net_Cooling"][:][:]
        data['metal_free_net_cooling'] = hf["/Metal_free/Net_Cooling"][:][:][:]
        data['metal_free_electron_density'] = hf["/Metal_free/Electron_density_over_n_h"][:]
        data['solar_electron_density'] = hf["/Solar/Electron_density_over_n_h"][:]
        data['He_massfrac_bins'] = hf["/Metal_free/Helium_mass_fraction_bins"][:]
    return data

def interpolate_cooling_table(data_file, density, temperature, H, He, metallicity):
    Mf_net_cool = data_file['metal_free_net_cooling'] #He bins (7), temp bins (176), density bins (41)
    Mf_ne = data_file['metal_free_electron_density'] #He bins (7), temp bins (176), density bins (41)
    M_net_cool = data_file['total_metals_net_cooling'] #temp bins (176), density bins (41)
    M_ne =  data_file['solar_electron_density'] #temp bins (176), density bins (41)
    
    density_bins = np.log10(data_file['H_bins'])
    temperature_bins = np.log10(data_file['Temp_bins'])
    
    metal_free_net_cool = np.zeros((41,176))
    metal_free_electron_density = np.zeros((41,176))
    x = data_file['He_massfrac_bins']
    He_frac = He/(He+H)
    for i in range(0,41):
        for j in range(0,176):
            y = Mf_net_cool[:,j,i]
            interpolation = interpolate.interp1d(x, y)
            metal_free_net_cool[i,j] = interpolation(He_frac)
            y = Mf_ne[:,j,i]
            interpolation = interpolate.interp1d(x, y)
            metal_free_electron_density[i,j] = interpolation(He_frac)

    f = interpolate.interp2d(temperature_bins, density_bins, metal_free_net_cool, kind='cubic')
    metal_free_Lambda = f(temperature, density)

    f = interpolate.interp2d(temperature_bins, density_bins, metal_free_electron_density, kind='cubic')
    metal_free_ne = f(temperature, density)

    f = interpolate.interp2d(density_bins, temperature_bins, M_net_cool, kind='cubic')
    metal_Lambda = f(density, temperature)

    f = interpolate.interp2d(density_bins, temperature_bins, M_ne, kind='cubic')
    metal_ne = f(density, temperature)

    net_cool = metal_free_Lambda
    net_cool += (metal_free_ne/metal_ne) * metal_Lambda * metallicity
    
    # Add Compton cooling off the CMB
    STEFAN = 7.5657e-15 # erg cm^-3 K^-4
    C_speed = 2.9979e10 # cm s^-1
    ELECTRONMASS = 9.10953e-28 # g
    THOMPSON = 6.6524587e-25 # cm^2
    BOLTZMANN = 1.3806e-16 # erg K^-1
    
    #t_cmb = 2.728 * (1.0 + redshift)
    #comp_add = -1. * (4.0 * STEFAN * THOMPSON * (t_cmb**4) / (ELECTRONMASS * C_speed))
    #comp_add *= BOLTZMANN  * (t_cmb - 10**temperature) * metal_free_ne
    #net_cool += comp_add
    
    return net_cool

def compute_net_cooling_normal_opt(redshift=0, density=0.1, temperature=1e5, H=0.752, He=0.248, metallicity=0.01):
    """This routine computes net cooling for a given metallicity (solar relative abundances).
    
    It calculates the net cooling rate (Lambda/n_H^2 [erg s^-1 cm^3]) given an array of redshifts, densities,
    temperatures, and abundances using the cooling tables prepared by Wiersma, Schaye, & Smith (2009).
    
    This routine uses the "total_metals" array such that only a metallicity is needed, and relative metal
    abundances are assumed to be solar relative.
    
    Parameters
    ----------
    redshift: float / array
    The redshift (for collisional ionization, enter a dummy value)
    
    density : float / array
    The hydrogen number density (n_H [cm^-3])
    
    temperature : float / array
    The temperature (T [K])
    
    H/He : float
    Abundance of element H/He.
    
    metallicity: float
    Metallicity (relative to solar, thus Z = 1 performs the calculation for solar abundances).
    Note that the solar metal mass fraction for these tables is defined as Z = 1 - X - Y = 0.0129
    
    Returns
    -------
    
    Lambda: float / array
    Net cooling rates, in (Lambda_net/n_H^2 [erg cm^3 s^-1])
    
    """
    density = np.log10(density)
    temperature = np.log10(temperature)
    
    table_name_z1, table_name_z2, z_range = look_for_table_name(redshift)
    data1 = read_wiersma_cooling_table(table_name_z1)
    net_cool1 = interpolate_cooling_table(data1,density, temperature, H, He, metallicity)
    data2 = read_wiersma_cooling_table(table_name_z2)
    net_cool2 = interpolate_cooling_table(data2,density, temperature, H, He, metallicity)
    
    net_cool = np.array([net_cool1[0],net_cool2[0]])
    interpolation = interpolate.interp1d(z_range, net_cool)
    cooling_rate = interpolation(redshift)
    
    return cooling_rate







