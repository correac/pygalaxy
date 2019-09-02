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
        dz = (z[select_file[0]]-redshift)/(z[select_file[0]]-z[select_file[0]+1])
    else:
        output2 = dir+file[select_file[0]-1]
        dz = (z[select_file[0]]-redshift)/(z[select_file[0]]-z[select_file[0]-1])
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
        data['Solar_mass_fractions'] = hf["/Header/Abundances/Solar_mass_fractions"][:]
    return data

def compute_net_cooling_normal_opt(redshift=0, density=0.1, temperature=1e5, H=0.752, He=0.248, metallicity=0.01)
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

    table_name_z1, table_name_z2, dz = look_for_table_name(redshift)
    data1 = read_wiersma_cooling_table(table_name_z1)
    hhe_cool1 = data1['metal_free_net_cooling'] #He bins (7), temp bins (176), density bins (41)
    hhe_ne1 = data1['metal_free_electron_density'] #He bins (7), temp bins (176), density bins (41)
    metal_cool1 = data1['total_metals_net_cooling'] #temp bins (176), density bins (41)
    solar_ne_nh1 =  data1['solar_electron_density'] #temp bins (176), density bins (41)

    data2 = read_wiersma_cooling_table(table_name_z2)
    hhe_cool2 = data2['metal_free_net_cooling']
    hhe_ne2 = data2['metal_free_electron_density']
    metal_cool2 = data2['total_metals_net_cooling']
    solar_ne_nh2 =  data2['solar_electron_density']

    tbl_hhecool = hhe_cool1 * dz + hhe_cool2 * (1. - dz)
    tbl_hhene = hhe_ne1 * dz + hhe_ne2 * (1. - dz)
    tbl_solar_ne_nh = solar_ne_nh1 * dz + solar_ne_nh2 * (1. - dz)
    tbl_metalcool = metal_cool1 * dz + metal_cool2 * (1. - dz)

    tbl_log_dens = data1['H_bins']
    tbl_log_temperature = data1['Temp_bins']

    metal_free_net_cool = np.zeros((41,176))
    metal_free_electron_density = np.zeros((41,176))
    x = data1['He_massfrac_bins']
    He_frac = He/(He+H)
    for i in range(0,41):
        for j in range(0,176):
            y = tbl_hhecool[:,j,i]
            interpolation = interpolate.interp1d(x, y)
            metal_free_net_cool[i,j] = interpolation(He_frac)
            y = tbl_hhene[:,j,i]
            interpolation = interpolate.interp1d(x, y)
            metal_free_electron_density[i,j] = interpolation(He_frac)

    f = interpolate.interp2d(tbl_log_temperature, tbl_log_dens, metal_free_net_cool, kind='cubic')
    metal_free_Lambda = f(temperature, density)

    f = interpolate.interp2d(tbl_log_temperature, tbl_log_dens, metal_free_electron_density, kind='cubic')
    metal_free_ne = f(temperature, density)

    f = interpolate.interp2d(tbl_log_temperature, tbl_log_dens, tbl_metalcool, kind='cubic')
    metal_Lambda = f(temperature, density)

    f = interpolate.interp2d(tbl_log_temperature, tbl_log_dens, tbl_solar_ne_nh, kind='cubic')
    metal_ne = f(temperature, density)

    net_cool = metal_free_Lambda
    net_cool += (metal_free_ne/metal_ne) * metal_Lambda * metallicity
    return net_cool


