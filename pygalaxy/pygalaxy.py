import numpy as np
from argumentparser import ArgumentParser
from hot_halo import M200accr, fhotacc, fhot
from galaxy_gas_accretion import cooling_radius, galaxy_gas_accretion
from plotter import output_plots

class pygal:

    def __init__(
        self,
        redshift: float,
        halo_mass: float,
        output: str,
        option: str,
    ):
        """
        Parameters
        ----------
        redshift: unyt.array.unyt_quantity
        halo mass: unyt.array.unyt_quantity
        directory: str
        option: str
        """

        self.halo_mass = halo_mass
        self.redshift = redshift
        self.output_path = output
        self.output_option = option

        self.halo_gas_accretion = np.log10(M200accr(self.halo_mass,self.redshift,commah_flag=False))
        self.cooling_radius = cooling_radius(self.halo_mass, self.redshift)
        self.galaxy_gas_accretion = np.log10(galaxy_gas_accretion(self.halo_mass, self.redshift, self.cooling_radius))
        self.hot_gas_accretion = fhotacc(self.halo_mass, self.redshift)
        self.halo_hot_gas = fhot(self.halo_mass, self.redshift)

    def output_data(self):

        filename = self.output_path + "Correa_etal2018_data_redshift_{0:.3f}".format(self.redshift) + ".txt"
        header = "# Rates of gas accretion onto haloes and galaxies, cooling radius, " \
                 "fraction of gas accreted onto haloes in hot mode, and fraction hot halo gas. \n"
        header += "# The derivation of data and routines follows the work of Correa et al. (2018a,b); MNRAS " \
                  "478, 1, p.255 and MNRAS 473, 1, p.538. \n"
        header += "# The data is output for redshift {0:.3f}".format(self.redshift) + ".\n"

        fout = open(filename, 'w')
        fout.write(header)

        fout.write("# Halo mass    -    Halo gas    -   Galaxy gas   - Cooling  - Hot mode - Hot halo gas" + '\n' + \
                   "#              -    Accretion   -   accretion    -  radius  - fraction -   fraction" + '\n' + \
                   "# log10[Msun]  - log10[Msun/yr] - log10[Msun/yr] - [/R200]  -          -  " + '\n')

        for i in range(len(self.halo_mass)):
            output = "   {0:0.3f}".format(self.halo_mass[i]) + "           {0:5.3f}".format(self.halo_gas_accretion[i])\
                     + "           {0:5.3f}".format(self.galaxy_gas_accretion[i]) + "       {0:.3f}".format(self.cooling_radius[i])\
                     + "      {0:.3f}".format(self.hot_gas_accretion[i]) + "      {0:.3f}".format(self.halo_hot_gas[i]) + '\n'
            fout.write(output)

        fout.close()  # Close file


if __name__ == "__main__":

    config_parameters = ArgumentParser()

    # Halo mass range [log10 Msun]:
    Mhalo = np.arange(11.0, 14.2, 0.2)

    # Load input options and calculate data
    pygalaxy = pygal(
        redshift=config_parameters.redshift,
        halo_mass=Mhalo,
        output=config_parameters.output_directory,
        option=config_parameters.output_option,
    )

    # Output .txt data file
    pygalaxy.output_data()

    # Output plots showing relations
    output_plots(pygalaxy)
