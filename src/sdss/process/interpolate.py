"""
Functionality to interpolate spectra in a common grid.
Functionality to do the computations in parallel
"""

import multiprocessing as mp
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from sdss.metadata import MetaData
from sdss.utils.managefiles import FileDirectory
from sdss.utils.parallel import to_numpy_array


class Interpolate(FileDirectory, MetaData):
    """
    Interpolate SDSS spectra in a common grid, follwing these steps

    1. Remove wavelengths where sky is problematic
    2. Remove large relative uncertainties (where std>flux)
    3. Deredenning spectrum
    4. Deredshift
    5. Interpolate

    NOTE: remove means to replace with a NaNs

    """

    def __init__(
        self,
        meta_data_df: pd.DataFrame,
        raw_data_dir: str,
        wave_parameters: dict,
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
        number_processes: number of jobs when processing a bulk
            of a spectra
        OUTPUT
            check how to document the constructor of a class
        """

        FileDirectory.__init__(self)
        MetaData.__init__(self)

        super().check_directory(raw_data_dir, exit_program=True)
        self.spectra_directory = raw_data_dir
        self.meta_data = meta_data_df
        self.grid = self.get_grid(wave_parameters)

        self.extinction = self.dust_model()

    def get_grid(self, wave_parameters: dict) -> np.array:
        """
        Computes the master grid for the interpolation of the spectra
        ARGUMENTS
            grid_parameters: dictionary with structure
                {
                    "number_waves": number fluxes in the grid,
                    "lower": lower bound in the grid,
                    "upper": upper bound in the grid
                }
        RETURN
            wave_grid: numpy array with the grid
        """

        number_waves = wave_parameters["number_waves"]
        lower = wave_parameters["lower"]
        upper = wave_parameters["upper"]

        grid = np.linspace(lower, upper, number_waves)

        return grid

    @staticmethod
    def OI_5577_interpolation(wave: np.array, spectrum: np.array) -> np.array:
        """
        Interpolate region of airglow radiation from the
        [OI]5577 region with the average flux of neighboring
        segments to that region

        INPUTS
        wave: wavelengths of the raw spectrum in observers
            frame
        spectrum: flux of spectrum in observers frame

        OUTPUT
        spectrum: flux with interpolation of the [OI]5577
            region caused by atmospheric airglow
        """

        OI_mask = np.logical_and(wave > 5565, wave < 5590)
        n_mask = OI_mask.sum()

        try:

            left_idx = np.argwhere(wave == wave[OI_mask][0])[0, 0]
            right_idx = np.argwhere(wave == wave[OI_mask][-1])[0, 0]

        except IndexError:

            return spectrum

        left_flux = spectrum[left_idx - n_mask : left_idx]
        right_flux = spectrum[right_idx : right_idx + n_mask]

        average_flux = (left_flux + right_flux) / 2.0

        spectrum[OI_mask] = average_flux

        return spectrum

    def interpolate(self, specobjid: int) -> tuple:
        """
        Interpolate a single spectrum

        INPUTS
        specobjid: specobj of a spectrum, the name of the file
            with raw data is f'{raw_data_directory}/{specobjid}.npy'

        OUTPUT
        spectrum, variance: interpolated spectrum and its variance
            over the common grid
        """

        spectrum_location = f"{self.spectra_directory}/{specobjid}.npy"

        spectrum = np.load(spectrum_location)

        wave = spectrum[0]
        flux = spectrum[1]
        ivar = spectrum[2]

        # remove [OI]5577 line
        flux = self.OI_5577_interpolation(wave=wave, spectrum=flux)
        # remove large uncertainties
        flux, variance = self.remove_large_uncertainties(flux, ivar)
        # correct for extinction
        ebv = self.meta_data.loc[specobjid, "ebv"]
        flux = self.dered_spectrum(flux, wave, ebv)

        # deredshift
        z = self.meta_data.loc[specobjid, "z"]
        wave = self.convert_to_rest_frame(wave, z)

        # interpolate to common grid
        flux = np.interp(self.grid, wave, flux, left=np.nan, right=np.nan)

        variance = np.interp(
            self.grid, wave, variance, left=np.nan, right=np.nan
        )

        return flux, variance

    def dered_spectrum(
        self, flux: np.array, wave: np.array, ebv: float
    ) -> np.array:

        """
        Apply extinction function to spectrum for deredening.

        INPUTS
        flux: spectrum fluxes in observer frame
        wave: wavelengths of spectrum in observer frame
        ebv: E(B-v) from from Schlegel, Finkbeiner & Davis (1998)

        OUTPUT
        dered_flux: fluxes corrected by extinction
        """

        # extinction array for this spectrum
        extinction = self.extinction(wave)
        extinction *= ebv

        dereded_flux = flux * 10 ** (extinction / 2.5)

        return dereded_flux

    def dust_model(self) -> object:
        """Extinction function to dered spectra"""
        # dust model for dereden_spectrum
        wave = np.array([2600, 2700, 4110, 4670, 5470, 6000, 12200, 26500])

        extinction = np.array(
            [6.591, 6.265, 4.315, 3.806, 3.055, 2.688, 0.829, 0.265]
        )

        extinction = interp1d(wave, extinction, kind="cubic")

        return extinction

    def remove_large_uncertainties(
        self, flux: np.array, ivar: np.array
    ) -> np.array:
        """
        Set to NaN fluxes with large variance in their measurements

        INPUT
        flux: raw spectrum
        ivar: inverse of variance for fluxe measurements of spectrum

        OUTPUTS

        flux, variance: masked flux with large noise in its measurements
            and variance of flux measurements
        """

        # Get variance of each flux
        ivar[ivar == 0] = np.nan
        variance = 1 / ivar
        variance[np.isnan(ivar)] = np.inf
        variance = np.nan_to_num(variance)

        # mask of: variance > flux [higly uncertain values]
        flux_no_nans = np.nan_to_num(flux)
        large_variance_mask = np.sqrt(variance) > flux_no_nans

        flux[large_variance_mask] = np.nan

        return flux, variance

    def convert_to_rest_frame(
        self,
        wave: np.array,
        z: float,
    ) -> np.array:

        """
        Convert from observer frame to rest frame according to
        redshift of galaxy

        INPUTS
        wave: wavelegths in observer frame
        z: redshift

        OUTPUTS
        wave: wavelengths in rest frame
        """

        rest_frame_factor = 1.0 / (1.0 + z)
        wave = wave * rest_frame_factor

        return wave


def shared_data(
    input_counter: mp.Value,
    input_meta_data: pd.DataFrame,
    input_grid_parameters: dict,
    input_raw_data_directory: str,
    shared_arrays_parameters: tuple,
) -> None:

    """
    Data to share with child processes during interpolation.

    INPUTS
    input_counter: value with lock to track each spectrum
    input_meta_data: data frame with meda data of all spectra
    input_grid_parameters:
        {
            "upper": upper bound in wavelength common grid,
            "lower": lower bound in wavelength common grid,
            "number_waves": number of elements in common grid
        }
    input_raw_data_directory: path to raw data
    share_arrays_parameters: contains data of shared arrays
        spectra: RawArray for spectra
        spectra_shape: (number_of_spectra, number_of_waves)
        variance: RawArray for variance of spectra's fluxes
        variance_shape: (number_of_spectra, number_of_waves)
        ids: RawArray to link the position of a spectrum with
            its specobjid in the spectra array
        ids_shape: (number_of_spectra, 2)
            column 0: spectra id, column 1: specobjid

    """

    global counter
    global meta_data
    global grid_parameters
    global raw_data_directory
    global spectra
    global variance_of_spectra
    global track_indexes
    global interpolator

    counter = input_counter
    meta_data = input_meta_data
    grid_parameters = input_grid_parameters
    raw_data_directory = input_raw_data_directory

    spectra = shared_arrays_parameters[0]
    spectra_shape = shared_arrays_parameters[1]
    spectra = to_numpy_array(spectra, spectra_shape)

    variance_of_spectra = shared_arrays_parameters[2]
    variance_shape = shared_arrays_parameters[3]
    variance_of_spectra = to_numpy_array(variance_of_spectra, variance_shape)

    track_indexes = shared_arrays_parameters[4]
    ids_shape = shared_arrays_parameters[5]
    track_indexes = to_numpy_array(track_indexes, ids_shape)

    interpolator = Interpolate(
        meta_data_df=meta_data,
        raw_data_dir=raw_data_directory,
        wave_parameters=grid_parameters,
    )


def worker_interpolation(specobjid: int) -> None:

    """
    Worker to do interpolation of spectra in parallel.
    The workflow per spectrum is:

    1. Remove wavelengths where sky is problematic
    2. Remove large relative uncertainties (where std>flux)
    3. Deredenning spectrum
    4. Deredshift
    5. Interpolate

    NOTE: remove means to replace with a NaNs

    INPUTS
    specobjid: unique identifier of spectrum in data frame
        containing meta data
    """

    spectrum, variance_of_spectrum = interpolator.interpolate(
        specobjid=specobjid
    )

    with counter.get_lock():

        counter_value = counter.value
        counter.value += 1

        print(f"[{counter_value}] Interpolate {specobjid}", end="\r")

    spectra[counter_value, :] = spectrum
    variance_of_spectra[counter_value, :] = variance_of_spectrum

    index_track = np.array([counter_value, specobjid], dtype=np.uint)

    track_indexes[counter_value, :] = index_track
