import ctypes
from functools import partial
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import os

####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

####################################################################
# from src.process.worker import init_worker
################################################################################
class DataProcess:
    def __init__(
        self,
        galaxies_frame: "pd.df",
        number_processes: "int",
        grid_parameters: "dict",
        data_directory: "str",
        output_directory: "str",
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

        self.data_directory = data_directory
        self.spectra_directory = f"{data_directory}/rest_frame"

        self._check_directory(output_directory)
        self.output_directory = output_directory

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
    def interpolate(self):
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory

        OUTPUT
        """

        def to_numpy_array(shared_array, array_shape):
            """Create a numpy array backed by a shared memory Array."""
            share_array = np.ctypeslib.as_array(shared_array)
            return share_array.reshape(array_shape)

        def init_worker(share_array, array_shape):

            global fluxes

            fluxes = to_numpy_array(share_array, array_shape)

        print(f"Interpolate spectra...")

        number_spectra = self.frame.shape[0]
        number_waves = self.grid.size
        fluxes_shape = (number_spectra, number_waves)

        shared_fluxes = RawArray(ctypes.c_float, number_spectra * number_waves)

        galaxy_indexes = self.frame.index.values

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_worker,
            initargs=(shared_fluxes, fluxes_shape),
        ) as pool:

            results = pool.map(self._interpolate, galaxy_indexes)

        fluxes = to_numpy_array(shared_fluxes, fluxes_shape)
        save_to = f"{self.output_directory}/fluxes_interp.npy"
        np.save(save_to, fluxes)

        return fluxes

    ###########################################################################
    def _interpolate(self, galaxy_index: "int"):

        spectrum_name = self._get_name(galaxy_index)
        print(f"Interpolate {spectrum_name}", end="\r")

        spectrum_location = f"{self.spectra_directory}/{spectrum_name}.npy"
        spectrum = np.load(spectrum_location)

        flux = np.interp(
            self.grid,
            spectrum[0],  # wave
            spectrum[1],  # flux
            left=np.nan,
            right=np.nan,
        )

        # fluxes is the shared array made global in worker initializer
        # galaxy index makes shure that location in array is the same
        # as in the data frame :)
        fluxes[galaxy_index, :] = flux[:]

    ###########################################################################
    def _get_name(self, galaxy_index: "int"):

        galaxy = self.frame.iloc[galaxy_index]

        plate = f"{galaxy['plate']:04}"
        mjd = f"{galaxy['mjd']}"
        fiberid = f"{galaxy['fiberid']:04}"

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        return spectrum_name

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

###############################################################################
