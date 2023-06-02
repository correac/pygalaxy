import numpy as np
import commah
from wiersma_cooling_table import compute_net_cooling_normal_opt
from scipy.optimize import fsolve

def fhot(M200,z):
    """
    Function that returns the fraction of halo gas that is hot,
    given the total halo mass and redshift. This function follows eqs. (4-10)
    of Correa et al. (2018), MNRAS, 473, 1, 538.
    Gas (within 0.15Rvir < r < Rvir) is defined to be hot if t_cool > t_dyn
    (cooling time is longer than dynamical time). The fraction is normalized by
    the universal baryon fraction. This function is a best-fitting expression
    of hot gas fractions in halos from the EAGLE simulations.
    
    Parameters
    ----------
    M200: float
    Logarithmic halo mass in units of M200 (mass enclosed within the virial radius,
    R200 defined as the radius at which the mean density is 200 times the critical).
    
    z: float / array
    Redshift.
    
    Returns
    -------
    fhot : float / array
    Fraction of halo gas that is hot.
    fhot = M_hot / (Omega_b/Omega_m) M200
    """
    x = M200-12.
    z_tilde = np.log10(1.+z)
    if z<=2.0:
        alpha = -0.79+0.31*z_tilde-0.96*z_tilde**2
        beta = 0.52-0.57*z_tilde+0.85*z_tilde**2
        gamma = -0.05
    else:
        alpha = -0.38-1.56*z_tilde+1.17*z_tilde**2
        beta = 0.12+0.94*z_tilde-0.55*z_tilde**2
        gamma = -0.05
    return 10**(alpha+beta*x+gamma*x**2) #fhot

def Mhot(M200,z,f_baryon=0.04825/0.307):
    """
    Function that returns the logarithmic halo gas mass that is hot, given
    the total halo mass and redshift. This function follows eqs. (4-10)
    of Correa et al. (2018), MNRAS, 473, 1, 538.
    Gas (within 0.15Rvir < r < Rvir) is defined to be hot if t_cool > t_dyn
    (cooling time is longer than dynamical time). The fraction is normalized by
    the universal baryon fraction. This function is a best-fitting expression
    of hot gas fractions in halos from the EAGLE simulations.
    
    Parameters
    ----------
    M200: float
    Logarithmic halo mass in units of M200 (mass enclosed within the virial radius,
    R200 defined as the radius at which the mean density is 200 times the critical).
    
    z: float / array
    Redshift.
    
    f_baryon : float (optional)
    Universal baryon fraction.
    
    Returns
    -------
    Mhot : float / array
    Logarithmic halo gas mass that is hot.
    """
    Mhot = np.log10(fhot(M200,z))+np.log10(f_baryon)+M200
    return Mhot

def fhotacc(M200,z):
    """
    Function that returns the fraction of gas accreted hot as a function of
    halo mass (M200) and redshift. This function follows eqs. (11-15) of Correa
    et al. (2018), MNRAS, 473, 1, 538.
    Gas is defined to be accreted hot if immediately after crossing the virial
    radius the gas temperature is larger than 10^5.5 K. This function is a b
    est-fitting expression of hot gas accretion in halos from the EAGLE simulations.
    
    Parameters
    ----------
    M200: float
    Logarithmic halo mass in units of M200 (mass enclosed within the virial radius,
    R200 defined as the radius at which the mean density is 200 times the critical).
    
    z: float / array
    Redshift.
    
    Returns
    -------
    fhotacc : float / array
    Fraction of gas accreted hot.
    """
    z_tilde = np.log10(1.+z)
    if z<=2:
        M_half=-0.15+0.22*z+0.07*z**2
        alpha = -1.86*10**(-1.26*z_tilde+1.29*z_tilde**2)
    if (z>2.0)&(z<=4.0):
        M_half=-0.25+0.53*z-0.07*z**2
        alpha = -0.46*10**(0.81*z_tilde-0.42*z_tilde**2)
    if z>4:
        M_half = 0.72+0.01*z
        alpha = -1.07
    M_half += 12.0
    x = M200-M_half
    x = 10**x
    fhotacc = 1./(1.+x**alpha)
    return fhotacc

def rhocrit(z,Omega_m=0.307,h=0.6777,cgs=True):
    """
    Function that returns the universal critical density given the redshift, and cosmological
    parameters Omega matter and hubble parameter.
    
    Parameters
    ----------
    z: float / array
    Redshift.
    
    Omega_m: float, optional
    Universal mass density including ordinary mass (baryons) and dark matter. Optional input.
    Default is 0.307 as estimated by Planck15.
    
    h : float, optional
    Hubble parameter.  Optional input. Default is 0.6777 as estimated by Planck15.
    
    cgs : bool, optional
    If 'True' (default), the critical density is in units of gr/cm^3. If 'False', the
    critical density is in units of Msun/Mpc^3.
    
    Returns
    -------
    rhocrit : float / array
    Critical density.
    """
    Omega_l = 1.0-Omega_m
    H0 = h * 100. #km/s/Mpc
    if cgs:
        H0 /= 3.086e19 #1/Mpc -> 1/km
        G = 6.67408e-8 #cm^3 g^-1 s^-2
        rho0 = 3*H0**2/(8.*np.pi*G) #gr cm^-3
    else:
        G = 4.30091e-3 #pc Msun^-1 (km/s)**2
        G *= 1e-3 #pc->Mpc (Mpc Msun^-1 (km/s)**2
        rho0 = 3*H0**2/(8.*np.pi*G) #Msun Mpc^-3
    return rho0 * (Omega_m*(1.0+z)**3+Omega_l)

def M200accr(M200,z,cosmo='planck15',f_baryon=0.04825/0.307,commah_flag=True):
    """
    Function that calls the accretion routine from the COMMAH
    package, which calculates the accretion rate of a halo at any
    redshift 'z', given the halo total mass 'Mi' at redshift z.
    
    Parameters
    ----------
    M200 : float
    Logarithmic halo mass at redshift zero. Note, the halo is in
    units of M200 (mass enclosed within the virial radius, R200,
    defined as the radius at which the mean density is 200 times
    the critical).
    
    z : float / array
    Redshift to solve acc_rate.
    
    cosmo : dict, optional.
    Dictionary of cosmological parameters, similar in format to:
    {'N_nu': 0,'Y_He': 0.24, 'h': 0.702, 'n': 0.963,'omega_M_0': 0.275,
    'omega_b_0': 0.0458,'omega_lambda_0': 0.725,'omega_n_0': 0.0,
    'sigma_8': 0.816, 't_0': 13.76, 'tau': 0.088,'z_reion': 10.6}.
    Default: planck15 (cosmological parameters of Planck15).
    
    Returns
    -------
    dMdt : float / array
    Accretion rate [Msol/yr] at redshift 'z'.
    """
    if commah_flag:
        commah_run = commah.run(cosmo, z, 10**M200, z)
        output = commah_run['dMdt'].flatten()
        output *= f_baryon
    else:
        x = M200-12.0
        if z<=4.0:
            a = 0.83+0.553*z-0.0523*z**2
            b = 1.436-0.149*z+0.007*z**2
            c = -0.134+0.099*z-0.023*z**2
            output = 10**(a+b*x+c*x**2)
        else:
            a = 3.287-0.401*z+0.045*z**2
            b = 1.016+0.003*z+0.002*z**2
            output = 10**(a+b*x)
    return output

def Tvir(M200,z):
    """
    Function that calculates the halo virial temperature given the
    halo total mass (M200) and redshift.
    
    Parameters
    ----------
    M200 : float
    Logarithmic halo mass at redshift zero. Note, the halo is in
    units of M200 (mass enclosed within the virial radius, R200,
    defined as the radius at which the mean density is 200 times
    the critical).
    
    z : float / array
    Redshift.
    
    Returns
    -------
    Tvir : float / array
    Virial temperature.
    """
    Msun = 1.98847e33 #Msun->gr
    yr_to_sec = 3.154e7 #yr->sec
    mu = 0.59
    k_b = 1.3806485e-23 #m^2 kg s^-2 K^-1
    k_b *= (1e2)**2*1e3 #m->cm, kg->g
    mp = 1.6726e-24 # proton mass g
    G = 6.67408e-8 #cm^3 g^-1 s^-2
    R200 = 10**M200 * Msun / (4. * np.pi * 200. * rhocrit(z) / 3.)
    R200 = R200**(1./3.) #Virial radius
    Vc2 = 10**M200 * Msun * G / R200 #Maximum circular velocity: Vc2 = Vc^2(R200.)
    T = mu * mp * Vc2 / (2. * k_b)
    #T = 10**5.3 * (10**(M200-12.0))**(2./3.) * (1.+z)
    return T

def Gamma_heat(M200,z,cosmo='planck15',f_baryon=0.04825/0.307):
    """
    Function that calculates the virial heating rate as defined in
    eq. (18) of Correa et al. (2018), MNRAS, 473, 1, 538.
    
    Parameters
    ----------
    M200 : float
    Logarithmic halo mass at redshift zero. Note, the halo is in
    units of M200 (mass enclosed within the virial radius, R200,
    defined as the radius at which the mean density is 200 times
    the critical).
    
    z : float / array
    Redshift to solve acc_rate.
    
    cosmo : dict, optional.
    Dictionary of cosmological parameters.
    Default: planck15 (cosmological parameters of Planck15).
    
    f_baryon : float, optional.
    Universal baryon fraction.
    
    Returns
    -------
    Gamma_heat : float / array
    Virial heating rate [erg/s] at redshift 'z'.
    """
    Msun = 1.98847e33 #Msun->gr
    yr_to_sec = 3.154e7 #yr->sec
    mu = 0.59
    k_b = 1.3806485e-23 #m^2 kg s^-2 K^-1
    mp = 1.6726e-27 # proton mass kg
    cte = k_b / (mu * mp) #(m/s)^2 K^-1
    cte *= (1e2)**2 #m->cm, (cm/s)^2 K^-1
    Gamma = (3./2.) * cte * Tvir(M200,z)  #(cm/s)^2
    Gamma *= M200accr(M200,z,cosmo,f_baryon) #(cm/s)^2 Msun/yr
    Gamma *= Msun/yr_to_sec #erg/s
    Gamma *= ((2./3.)*fhot(M200,z) + fhotacc(M200,z))
    return Gamma

def Lambda(z, rho, T, metallicity, Xh=0.752, Xhe = 0.248):
    """
    This routine computes net cooling for a given redshift, density,
    temperature & metallicity (solar relative abundances) using the cooling
    tables prepared by Wiersma, Schaye, & Smith (2009).
    
    Parameters
    ----------
    redshift: float / array
    
    density : float / array
    The hydrogen number density (n_H [cm^-3])
    
    temperature : float / array
    The temperature (T [K])
    
    metallicity: float
    Metallicity (relative to solar, thus Z = 1 performs the calculation
    for solar abundances). Note that the solar metal mass fraction for
    these tables is defined as Z = 1 - X - Y = 0.0129
    
    XH/XHe : float, optional.
    Abundance of element H/He.
    
    Returns
    -------
    Lambda: float / array
    Net cooling rates, in (Lambda_net [erg cm^3 s^-1])
    """
    mu = 0.59
    mp = 1.6726e-24     # proton mass: gr
    n_gas = rho/(mu * mp) # gas number density [cm^-3]
    n_H = n_gas*Xh # hydrogen number density
    net_cool = compute_net_cooling_normal_opt(z, n_H, T, Xh, Xhe, metallicity) #return Lambda_net/n_H^2 [erg cm^3 s^-1]#
    net_cool *= n_H**2 #[erg cm^-3 s^-1]
    return net_cool

def Gamma_cool(M200,z,Thot=1e6,rho_hot=1.0,Zhot=0.1):
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
    rho_hot = 10**0.6*rhocrit(z) #density of hot gas (cgs=True: gr/cm^3)
    Gamma = 10**Mhot(M200,z) * Msun #units gr
    Gamma *= Lambda(z, rho_hot, Thot, Zhot) #units gr erg cm^-3 s^-1
    Gamma /= rho_hot #units erg/s
    return Gamma

def Mcrit(x,z):
    f = Gamma_heat(x,z)-Gamma_cool(x,z)
    return f

def calculate_M_critical(z):
    Mcritz = []
    initial_quess = 12.0
    for zi in z:
        Mi = fsolve(Mcrit, initial_quess, args=(zi))
        Mcritz = np.append(Mcritz,Mi)
        initial_guess = Mi
    
    return Mcritz
