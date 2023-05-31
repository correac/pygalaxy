import numpy as np
from argumentparser import ArgumentParser
from hot_halo import M200accr, fhotacc, fhot
from galaxy_gas_accretion import cooling_radius, galaxy_gas_accretion

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

        self.halo_gas_accretion = M200accr(self.halo_mass,self.redshift,commah_flag=False)
        self.cooling_radius = cooling_radius(self.halo_mass, self.redshift)
        self.galaxy_gas_accretion = galaxy_gas_accretion(self.halo_mass, self.redshift, self.cooling_radius)
        self.hot_gas_accretion = fhotacc(self.halo_mass, self.redshift)
        self.halo_hot_gas = fhot(self.halo_mass, self.redshift)

    


if __name__ == "__main__":

    config_parameters = ArgumentParser()

    # Halo mass range [log10 Msun]:
    Mhalo = np.arange(11, 12.2, 0.2)

    # Load input options and calculate data
    pygalaxy = pygal(
        redshift=config_parameters.redshift,
        halo_mass=Mhalo,
        output=config_parameters.output_directory,
        option=config_parameters.output_option,
    )

    # pygalaxy.make_data()
    #
    # pygalaxy.make_plots()
    #
    # if pygalaxy.output_option == 'yes':
    #     pygalaxy.calculate_mcrit()