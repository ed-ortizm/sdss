import multiprocessing as mp
import os
import sys
import warnings

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from src.superclasses import FileDirectory
from src.superclasses import MetaData
###############################################################################
def init_worker(
    input_counter: "mp.Value", input_df: "pandas dataframe"
) -> "None":
    """
    Initialize worker
    PARAMETERS
        counter: counts the number of the child process
        input_df: pandas dataframe to be accessible to each child
    """
    global counter
    global files_df

    counter = input_counter
    files_df = input_df

###############################################################################
class RawData(FileDirectory, MetaData):
    """Get wave, flux and ivar from sdss dr16 spectra"""

    def __init__(
        self,
        data_directory: "str",
        output_directory: "str",
        number_processes: "int",
    ):
        """
        PARAMETERS

            data_directory : sdss raw data's directory
            output_directory : save here .npy files
            number_processes : number of processes to use with mp.Pool

        OUTPUT
            RawData object
        """

        self.number_processes = number_processes

        super().check_directory(data_directory, exit=True)
        self.data_directory = data_directory

        super().check_directory(output_directory)
        self.output_directory = output_directory


    ###########################################################################
    def remove_fits_files(self, files_df: "pandas_data_frame")->"None":
        """
        Remove unwanted files, e.g. non galaxies fro sample

        PARAMETER
            files_df: data frame with necessary data to delete
                unwanted files
        """
        print(f"Remove files...")

        # I use specobjid as the index in the data frame
        files_indexes = files_df.index.values

        counter = mp.Value("i", 0)

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_worker,
            initargs=(counter, files_df),
        ) as pool:

            pool.map(self._remove_fits_file, files_indexes)

        print(f"Remove files finish...")
    ###########################################################################
    def _remove_fits_file(self, file_index: "int")->"None":
        """
        Remove file
        PARAMETERS
            file_index: specobjid of the object
        """

        file_row = files_df.loc[file_index]

        [file_directory, spectrum_name] = self._get_file_location(file_row)

        file_location = f"{file_directory}/{spectrum_name}.fits"

        with counter.get_lock():
            counter.value += 1
            print(f"[{counter.value}] Remove {spectrum_name}", end="\r")

        super().remove_file(file_location)


    ###########################################################################
    def save_raw_data(self, files_df: "pandas dataframe") -> "None":
        """
        Save data frame with all meta dat

        PARAMETERS
            files_df: data frame with meta data of galaxies
        """

        print(f"Save wave, flux, ivar!")

        files_indexes = files_df.index.values

        counter = mp.Value("i", 0)

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_worker,
            initargs=(counter, files_df),
        ) as pool:

            results = pool.map(self._get_data, files_indexes)

        number_fail = sum(results)
        print(f"Fail with {number_fail} files")

    ###########################################################################
    def _get_data(self, file_index: "int") -> "int":
        """
        Save data from spetrum corresponding to file_index in
        files_df

        PARAMETERS
            file_index: specobjid of a galaxy in that galaxy that the
                frame passed to the constructor of the class

        OUTPUT
            0 for successful operation, 1 otherwise
        """

        file_row = files_df.loc[file_index]

        [
            sas_directory,
            spectrum_name
        ] = super().get_file_location_sas(file_row)

        # [file_directory, spectrum_name] = self._get_file_location(file_row)

        file_location = (
            f"{self.data_directory}/{sas_directory}/{spectrum_name}.fits"
        )

        if not super().file_exists(file_location, exit=False):
            print(file_location)

            return 1

        with counter.get_lock():
            counter.value += 1
            print(f"[{counter.value}] Get {spectrum_name}", end="\r")

        result = self._get_save_wave_flux_ivar(
            file_index, file_location, spectrum_name
        )

        return result

    ###########################################################################
    def _get_save_wave_flux_ivar(
        self,
        file_index: "int",
        file_location: "str",
        spectrum_name: "str",
        # number_attempts: "int" = 10
    ) -> "int":

        """
        Get and save wav, flux and ivar
        PARAMETER

            file_index: specobjid of the galaxy in the input data frame
                passed to the constructor
            file_location: location of the fits file including extension
            spectrum_name: name of spectrum
                spec-{plate}-{mjd}-{fiberid}
            number_attempts: [no implement yet]
                how many times to try if warning

        OUTPUT
            0 for successful operation, 1 otherwise
        """

        save_to = f"{self.output_directory}/{spectrum_name}.npy"

        if super().file_exists(save_to, exit=False):
            print(f"Data of {spectrum_name} already saved!", end="\r")
            return 0

        warnings.filterwarnings(action="error")

        try:

            hdul = pyfits.open(file_location, memmap=False)

            wave = np.array(10.0 ** (hdul[1].data["loglam"]), dtype=float)
            flux = hdul[1].data["flux"]
            ivar = hdul[1].data["ivar"]
            specobjid = int(hdul[2].data["specobjid"].item())

            assert specobjid == file_index, "specobjid do not match"

            hdul.close()

            array_to_save = np.vstack((wave, flux, ivar))

            np.save(save_to, array_to_save)

            return 0

        except Exception as e:

            print(f"Problem with {spectrum_name}")
            print(e)

            return 1

###############################################################################
