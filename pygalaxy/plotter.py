import matplotlib.pyplot as plt
from pylab import *
import numpy as np

def output_plots(pygalaxy):

    #################
    # Plot parameters
    params = {
        "font.size": 14,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (8, 6),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.98,
        "figure.subplot.bottom": 0.10,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 2,
        "lines.linewidth": 1.5,
    }
    plt.rcParams.update(params)

    plt.figure()
    ax = plt.subplot(2, 2, 1)
    plt.grid(linestyle='-', linewidth=0.3)

    plt.plot(pygalaxy.halo_mass, pygalaxy.halo_gas_accretion, '-', color='tab:blue')

    plt.ylabel(r"$\log_{10}~\dot{M}_{\mathrm{gas,halo}}$ [M$_{\odot}$/yr]")
    plt.xlabel(r"$\log_{10}~M_{\mathrm{200c}}$ [M$_{\odot}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    #####
    ax = plt.subplot(2, 2, 2)
    plt.grid(linestyle='-', linewidth=0.3)

    plt.plot(pygalaxy.halo_mass, pygalaxy.halo_gas_accretion, '--', color='tab:blue')
    plt.plot(pygalaxy.halo_mass, pygalaxy.galaxy_gas_accretion, '-', color='tab:blue')

    plt.ylabel(r"$\log_{10}~\dot{M}_{\mathrm{gas,galaxy}}$ [M$_{\odot}$/yr]")
    plt.xlabel(r"$\log_{10}~M_{\mathrm{200c}}$ [M$_{\odot}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    #####
    ax = plt.subplot(2, 2, 3)
    plt.grid(linestyle='-', linewidth=0.3)

    plt.plot(pygalaxy.halo_mass, pygalaxy.cooling_radius, '-', color='tab:blue')

    plt.ylabel(r"$r_{\mathrm{cool}}/R_{\mathrm{200c}}$")
    plt.xlabel(r"$\log_{10}~M_{\mathrm{200c}}$ [M$_{\odot}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    #####
    ax = plt.subplot(2, 2, 4)
    plt.grid(linestyle='-', linewidth=0.3)

    plt.plot(pygalaxy.halo_mass, pygalaxy.hot_gas_accretion, '-', color='tab:blue')

    plt.ylabel(r"$f_{\mathrm{accr,hot}}$")
    plt.xlabel(r"$\log_{10}~M_{\mathrm{200c}}$ [M$_{\odot}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    filename = pygalaxy.output_path + "Correa_etal2018_data_redshift_{0:.3f}".format(pygalaxy.redshift) + ".png"
    plt.savefig(filename, dpi=300)
    plt.close()
