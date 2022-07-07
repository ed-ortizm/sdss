import ctypes
from functools import partial
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import os
import warnings

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from sdss.utils.managefiles import FileDirectory
from sdss.metadata import MetaData

###############################################################################
def to_numpy_array(input_shared_array, array_shape):
    """Create a numpy array backed by a shared memory Array."""

    share_array = np.ctypeslib.as_array(input_shared_array)

    return share_array.reshape(array_shape)


def init_process_worker(
    input_counter: mp.Value,
    input_df: pd.DataFrame,
    input_share_array: RawArray,
    array_shape: tuple,
    input_share_track_indexes: RawArray,
    input_share_track_indexes_shape: tuple,
) -> None:
    """
    Initialize worker to get sample relevant for the science
    PARAMETERS
        counter: counts the number of the child process
        input_df: pandas dataframe to be accessible to each child
        input_share_array:
        array_shape:
        input_share_track_indexes:
    """
    global counter
    global spectra_df
    global fluxes
    global track_indexes

    counter = input_counter
    spectra_df = input_df

    fluxes = to_numpy_array(input_share_array, array_shape)

    track_indexes = to_numpy_array(
        input_share_track_indexes, input_share_track_indexes_shape
    )


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
    def interpolate(self, spectra_df: pd.DataFrame) -> np.array:
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory
        OUTPUT
            fluxes: contains interpolate fluxes
        """

        number_spectra = spectra_df.shape[0]
        print(f"Interpolate {number_spectra} spectra")
        number_waves = self.grid.size

        counter = mp.Value("i", 0)

        fluxes = RawArray(ctypes.c_double, number_spectra * number_waves)
        fluxes_shape = (number_spectra, number_waves)

        # array with counter, specobjid and 0 or 1 if processing is
        # successful columns
        track_indexes = RawArray(ctypes.c_uint64, number_spectra * 3)
        track_indexes_shape = (number_spectra, 3)

        spectra_indexes = spectra_df.index.to_numpy()

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_process_worker,
            initargs=(
                counter,
                spectra_df,
                fluxes,
                fluxes_shape,
                track_indexes,
                track_indexes_shape,
            ),
        ) as pool:

            results = pool.map(self._interpolate, spectra_indexes)

        track_indexes = to_numpy_array(track_indexes, track_indexes_shape)
        success_interpolation_mask = ~track_indexes[:, 2].astype(np.bool)
        track_indexes = track_indexes[success_interpolation_mask]
        save_to = f"{self.output_directory}/indexes_interpolate.npy"
        np.save(save_to, track_indexes)

        fluxes = to_numpy_array(fluxes, fluxes_shape)

        fluxes = fluxes[success_interpolation_mask]

        save_to = f"{self.output_directory}/interpolate.npy"
        np.save(save_to, fluxes)

        return fluxes

    ###########################################################################
    def _interpolate(self, spectrum_index: "int") -> "None":
        """
        Interpolate a single spectrum
        PARAMETERS
            spectrum_index: specobj of a spectrum in spectra_df
        """

        spectrum_location = f"{self.spectra_directory}/{spectrum_index}.npy"

        try:

            spectrum = np.load(spectrum_location)

            wave = spectrum[0]
            flux = spectrum[1]
            # remove sky emission
            flux[np.bitwise_and(wave > 5565, wave < 5590)] = np.nan

            z = spectra_df.loc[spectrum_index, "z"]
            wave = self._convert_to_rest_frame(wave, z)

            flux = np.interp(self.grid, wave, flux, left=np.nan, right=np.nan)

            with counter.get_lock():

                fluxes[counter.value, :] = flux[:]

                index_track = np.array(
                    [counter.value, spectrum_index, 0], dtype=np.uint
                )
                track_indexes[counter.value, :] = index_track[:]

                print(
                    f"[{counter.value}] Interpolate {spectrum_index}", end="\r"
                )

                counter.value += 1

            return 0

        except Exception as e:

            with counter.get_lock():
                # ?
                dummy_flux = np.empty(fluxes.shape[1]) * np.nan

                fluxes[counter.value, :] = dummy_flux[:]

                index_track = np.array(
                    [counter.value, spectrum_index, 1], dtype=np.uint
                )
                track_indexes[counter.value, :] = index_track[:]

                print(
                    f"[{counter.value}] Fail to interpolate {spectrum_index}",
                    end="\r",
                )

                counter.value += 1

                print(e, index_track)

            return 1

    ###########################################################################
    def _convert_to_rest_frame(self, wave: "np.array", z: "float"):

        rest_frame_factor = 1.0 / (1.0 + z)
        wave = wave * rest_frame_factor

        return wave

    ###########################################################################
    def replace_missing_fluxes_and_normalize_by_median(
        self, spectra: "np.array"
    ) -> "np.array":
        """ """

        missing_values_mask = ~np.isfinite(spectra)

        nan_median = np.nanmedian(spectra, axis=1)
        null_median_mask = nan_median == 0.0
        nan_median[null_median_mask] = 1.0

        spectra *= 1 / nan_median.reshape((-1, 1))

        spectra[missing_values_mask] = 1.0

        # zero for spectra where nanmedian is null
        spectra[null_median_mask, :] = 0

        return spectra

    ###########################################################################
    def drop_indefinite_values(
        self, spectra: np.array, drop: float = 0.1
    ) -> list:
        """
        Drops indefinite values per wavelength according to
        the drop fraction specified by the drop parameter
        PARAMETERS
            spectra: contains interpolated fluxes in a common grid
            drop: fraction of indefinite values to discard at a given
                wavelength, e.g. drop = 0.1 will discard a wavelength
                where more than 10% of fluxes are indefinite values
        OUTPUTS
            [spectra, wave]:
                spectra: fluxes without indefinite values
                wave: common grid to all spectra update according
                    to wavelengths dropped
        """
        print(f"Shape before discard indefinite values: {spectra.shape}")

        number_indefinite_fluxes = np.count_nonzero(
            ~np.isfinite(spectra), axis=0
        )

        number_fluxes = spectra.shape[0]
        keep_fluxes_mask = number_indefinite_fluxes < number_fluxes * drop

        spectra = spectra[:, keep_fluxes_mask]

        print(f"Shape after discard indefinite values: {spectra.shape}")

        # update wavelength grid
        wave = self.grid[keep_fluxes_mask]

        return [spectra, wave]