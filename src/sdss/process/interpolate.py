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


class DataProcess(FileDirectory, MetaData):
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
        meta_data: pd.DataFrame,
        raw_data_directory: str,
        grid_parameters: dict,
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

        FileDirectory.__init__(self)
        MetaData.__init__(self)
        
        super().check_directory(raw_data_directory, exit_program=True)
        self.spectra_directory = raw_data_directory
        self.meta_data = meta_data
        self.grid = self.get_grid(grid_parameters)
        
        self.extinction = self.dust_model()


    def get_grid(self, grid_parameters: dict) -> np.array:
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

    def interpolate(self, specobjid: int) -> tuple(np.array, np.array):
        """
        Interpolate a single spectrum
        PARAMETERS
        spectrum_index: specobj of a spectrum, the name of the file
            with raw data is f'{raw_data_directory}/{specobjid}.npy'
        """

        spectrum_location = f"{self.spectra_directory}/{specobjid}.npy"

        spectrum = np.load(spectrum_location)

        wave = spectrum[0]
        flux = spectrum[1]
        ivar = spectrum[2]

        # remove sky emission
        flux[np.bitwise_and(wave > 5565, wave < 5590)] = np.nan
        # remove large uncertainties
        flux, variance = self.remove_large_uncertanties(flux, ivar)
        # correct for extinction
        ebv = self.meta_data.loc[specobjid, "ebv"]
        flux = self.dered_spectrumb(flux, wave, ebv)

        # deredshift
        z = self.meta_data.loc[specobjid, "z"]
        wave = self.convert_to_rest_frame(wave, z)

        # interpolate to common grid
        flux = np.interp(self.grid, wave, flux, left=np.nan, right=np.nan)
        
        variance = np.interp(
            self.grid, wave, variance, left=np.nan, right=np.nan
        )
        
        return flux, variance

    def dered_spectrum(self,
        flux: np.array,
        wave: np.array,
        ebv: float
    )-> np.array:
        
        # extinction array for this spectrum
        extinction = self.extinction(wave)
        extinction *= ebv

        dereded_flux = flux * 10 ** (extinction/2.5)

        return dereded_flux
        
    def dust_model(self) -> object:
        # dust model for dereden_spectrum
        wave = np.array(
            [ 2600,  2700,  4110,  4670,  5470,  6000, 12200, 26500]
        )
        
        extinction = np.array(
            [ 6.591,  6.265,  4.315,  3.806,  3.055,  2.688,  0.829,  0.265]
        )
        
        extinction = interp1d(wave, extinction, kind="cubic")
        
        return extinction
        
    def remove_large_uncertainties(self,
        flux: np.array, ivar: np.array
    )-> tuple(np.array, np.array):
        
        # Get variance of each flux
        ivar[ivar == 0] = np.nan
        variance = 1/ivar
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

        rest_frame_factor = 1.0 / (1.0 + z)
        wave = wave * rest_frame_factor

        return wave


def worker_interpolation(specobjid: int):

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
