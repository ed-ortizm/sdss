import ctypes
from functools import partial
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import os

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from src.superclasses import FileDirectory
from src.superclasses import MetaData
###############################################################################
def to_numpy_array(input_shared_array, array_shape):
    """Create a numpy array backed by a shared memory Array."""

    share_array = np.ctypeslib.as_array(input_shared_array)

    return share_array.reshape(array_shape)

def init_process_worker(
    input_counter: "mp.Value",
    # input_df: "pandas dataframe",
    input_share_array: "c types share array",
    array_shape: "tuple",
) -> "None":
    """
    Initialize worker to get sample relevant for the science
    PARAMETERS
        counter: counts the number of the child process
        input_df: pandas dataframe to be accessible to each child
    """
    global counter
    # global spectra_df
    global share_array

    counter = input_counter
    # spectra_df = input_df
    share_array = to_numpy_array(input_share_array, array_shape)

###############################################################################
class DataProcess(FileDirectory, MetaData):
    """Process sdss spectra"""
    def __init__(
        self,
        raw_data_directory: "str",
        output_directory: "str",
        grid_parameters: "dict",
        number_processes: "int",
    ):
        """
        Class to process  spectra

        PARAMETERS
            data_directory: location of raw spectra
            output_directory:
            grid_parameters: dictionary with structure
                {
                    "number_waves": "number fluxes in the grid",
                    "lower": "lower bound in the grid",
                    "upper": "upper bound in the grid"
                }
            number_processes: number of jobs when processing a bulk of a spectra
        OUTPUT
            check how to document the constructor of a class
        """

        super().check_directory(raw_data_directory, exit=True)
        self.spectra_directory = raw_data_directory

        super().check_directory(output_directory, exit=False)
        self.output_directory = output_directory

        self.grid = self._get_grid(grid_parameters)

        self.number_processes = number_processes
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
    def interpolate(self, spectra_df: "pandas dataframe")->"np.array":
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory

        OUTPUT
            share_array: contains interpolate fluxes
        """


        number_spectra = spectra_df.shape[0]
        print(f"Interpolate {number_spectra} spectra")
        number_waves = self.grid.size

        counter = mp.Value("i", 0)

        fluxes = RawArray(ctypes.c_float, number_spectra * number_waves)
        fluxes_shape = (number_spectra, number_waves)


        spectra_indexes = spectra_df.index

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_process_worker,
            initargs=(counter, fluxes, fluxes_shape),
        ) as pool:

            results = pool.map(self._interpolate, spectra_indexes)

        fluxes = to_numpy_array(fluxes, fluxes_shape)

        save_to = f"{self.output_directory}/interpolate.npy"
        np.save(save_to, fluxes)

        return fluxes

    ###########################################################################
    def _interpolate(self, spectrum_index: "int")->"None":
        """
        Interpolate a single spectrum

        PARAMETERS
            spectrum_index: specobj of a spectrum in spectra_df

        """

        # spectrum_row_df = spectra_df.loc[spectrum_index]

        spectrum_location = f"{self.spectra_directory}/{spectrum_index}.npy"
        spectrum = np.load(spectrum_location)

        with counter.get_lock():

            flux = np.interp(
                self.grid,
                spectrum[0],  # wave
                spectrum[1],  # flux
                left=np.nan,
                right=np.nan,
            )

            #?
            share_array[counter.value, :] = flux[:]

            counter.value += 1
            print(f"[{counter.value}] Interpolate {spectrum_index}", end="\r")

    ###########################################################################
    def normalize(self, spectra: "np.array"):
        """Spectra has no missing values"""

        normalization_array = np.median(spectra, axis=1)
        normalization_array = normalization_array.reshape((-1, 1))

        spectra[:, :] *= 1 / normalization_array

        return spectra

    ###########################################################################
    def replace_missing_flux(self, spectra: "array", method: "str" = "median"):
        """"""
        #######################################################################
        if method == "median":

            mask_replacement = ~np.isfinite(spectra)

            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmedian(spectra[idx, :])

        elif method == "mean":

            mask_replacement = ~np.isfinite(spectra)
            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmean(spectra[idx, :])

        return spectra

    ###########################################################################
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
