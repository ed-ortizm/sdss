import ctypes
from functools import partial
import multiprocessing as mp
import os

####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

################################################################################
class DataProcess:
    def __init__(
        self,
        galaxies_frame: "pd.df",
        number_processes: "int",
        grid_parameters: "dict",
    ):
        """
        Class to process rest frame spectra
        PARAMETERS
            galaxy_frame: meta data of sdss galaxies, such as name,
                signal to noise ratio and z
            number_processes: number of jobs when processing a bulk of a spectra
        OUTPUT
            check how to document the constructor of a class
        """

        self.frame = galaxies_frame
        self.number_processes = number_processes
        self.grid = self._get_grid(grid_parameters)

        self.fluxes = None

    ###########################################################################
    def interpolate_parallel(
        self, data_directory: "str", output_directory: "str"
    ):

        ######################################################################
        def to_numpy_array(shared_array, shape):
            """Create a numpy array backed by a shared memory Array."""
            fluxes = np.ctypeslib.as_array(shared_array)
            return fluxes.reshape(shape)

        def init_worker(shared_array, shape):
            """
            Initialize worker for processing:
            Create the numpy array from the shared memory Array for each process in the pool.
            """
            global fluxes
            fluxes = to_numpy_array(shared_array, shape)

        ######################################################################
        spectra_directory = f"{data_directory}/rest_frame"
        self._check_directory(spectra_directory)

        worker_function = partial(
            self._interpolate_parallel, spectra_directory=spectra_directory
        )
        ######################################################################
        print(f"Interpolate parallel...")

        number_spectra = self.frame.shape[0]
        number_fluxes = self.grid.size
        shape = (number_spectra, number_fluxes)

        shared_array = mp.Array(
            ctypes.c_double, number_spectra * number_fluxes, lock=False
        )

        fluxes = to_numpy_array(shared_array, shape)

        galaxy_indexes = self.frame.index

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_worker,
            initargs=(shared_array, shape),
        ) as pool:

            pool.map(worker_function, galaxy_indexes)

        self._check_directory(output_directory)
        self.frame.to_csv(f"{output_directory}/meta_data.csv", index=False)
        np.save(f"{output_directory}/fluxes_interp.npy", fluxes)

        return fluxes

    ###########################################################################
    def _interpolate_parallel(self, galaxy_index, spectra_directory):

        galaxy_name = self.frame.name.iloc[galaxy_index]
        print(f"Interpolate {galaxy_name}", end="\r")

        spectrum = np.load(f"{spectra_directory}/{galaxy_name}.npy")

        flux = np.interp(
            self.grid,
            spectrum[0],  # wave
            spectrum[1],  # flux
            left=np.nan,
            right=np.nan,
        )

        fluxes[galaxy_index, :] = flux[:]

    ###########################################################################
    def _check_directory(self, directory: "str", exit: "bool" = False):
        """
        Check if a directory exists, if not it creates it or
        exits depending on the value of exit
        """

        if not os.path.exists(directory):

            if exit:
                print(f"Directory {diretory} NOT FOUND")
                print("Code cannot execute")
                sys.exit()

            os.makedirs(directory)

    ###########################################################################
    def _get_grid(self, grid_parameters: "dict") -> "np.array":
        """
        Computes the master grid for the interpolation of the spectra

        ARGUMENTS

        grid_parameters: dictionary with structure
        {
        "number_waves": "number fluxes in the grid",
        "lower": "lower bound in the grid",
        "upper": "upper bound in the grid"
        }
        RETURN
        wave_grid: numpy array with the grid
        """
        number_waves = int(grid_parameters["number_waves"])
        lower = float(grid_parameters["lower"])
        upper = float(grid_parameters["upper"])

        grid = np.linspace(lower, upper, number_waves)

        return grid

    ###########################################################################
    def interpolate(self, data_directory: "str", output_directory: "str"):
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory

        PARAMETERS
            data_directory:
            output_directory:

        OUTPUT
            mp.pool list with integers telling whether the process was successful or not.
            Interpolate spectra will be saved in output directory
        """
        print(f"Interpolate spectra...")

        number_spectra = self.frame.shape[0]
        fluxes = np.empty((number_spectra, self.grid.size))

        galaxy_names = self.frame.name

        spectra_directory = f"{data_directory}/rest_frame"
        self._check_directory(spectra_directory)

        for idx, galaxy_name in enumerate(galaxy_names):

            flux = self._interpolate(galaxy_name, spectra_directory)
            fluxes[idx, :] = flux[:]

        self._check_directory(output_directory)
        self.frame.to_csv(f"{output_directory}/meta_data.csv", index=False)
        np.save(f"{output_directory}/fluxes_interp.npy", fluxes)

        return fluxes

    ###########################################################################
    def _interpolate(self, galaxy_name, spectra_directory):

        print(f"Interpolate {galaxy_name}", end="\r")

        spectrum = np.load(f"{spectra_directory}/{galaxy_name}.npy")

        flux = np.interp(
            self.grid,
            spectrum[0],  # wave
            spectrum[1],  # flux
            left=np.nan,
            right=np.nan,
        )

        return flux

    ###########################################################################
    def normalize(self, spectra: "np.array"):
        """Spectra has no missing values"""

        normalization_array = np.median(spectra, axis=1)
        normalization_array = normalization_array.reshape((-1, 1))

        spectra[:, :] *= 1 / normalization_array

        return spectra

    ############################################################################
    def replace_missing_flux(self, spectra: "array", method: "str" = "median"):
        """"""
        ########################################################################
        if method == "median":

            mask_replacement = ~np.isfinite(spectra)

            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmedian(spectra[idx, :])

        elif method == "mean":

            mask_replacement = ~np.isfinite(spectra)
            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmean(spectra[idx, :])

        return spectra

    ############################################################################
    def drop_indefinite_values(self, spectra: "np.array", drop: "float" = 0.1):
        """
        spectra: train
        discard_fraction:'float'=0.1
        """
        ########################################################################
        print(f"spectra shape before keep_spec_mask: {spectra.shape}")

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f"Indefinite vals in the input array: {np.sum(n_indef)}")

        keep_flux_mask = n_indef < spectra.shape[0] * drop

        spectra = spectra[:, keep_flux_mask]
        print(f"spectra shape after keep_spec_mask: {spectra.shape}")

        wave = self.grid[keep_flux_mask]

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f"Indefinite vals in the NEW array: {np.sum(n_indef)}")

        return spectra, wave

    ###########################################################################


###############################################################################
