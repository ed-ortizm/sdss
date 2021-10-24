import multiprocessing as mp
import os
import sys
import warnings

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from src.superclasses import FileDirectory

###############################################################################
# def init_get_data_worker(
#     input_counter: "mp.Value", input_df: "pandas dataframe"
# ) -> "None":
#     """
#     Initialize worker for download
#     PARAMETERS
#         counter: counts the number of the child process
#         input_df: pandas dataframe to be accessible to each child
#     """
#     global counter
#     global galaxies_df
#
#     counter = input_counter
#     galaxies_df = input_df
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
class GetRawData(FileDirectory):
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
            GetRawData object
        """

        self.number_processes = number_processes

        super().check_directory(data_directory, exit=True)
        self.data_directory = data_directory

        super().check_directory(output_directory)
        self.data_output_directory = output_directory

        super().check_directory(f"{output_directory}/rest_frame")
        self.rest_frame_directory = f"{output_directory}/rest_frame"

    ###########################################################################
    def remove_fits_files(self, files_df: "pandas_data_frame")->"None":
        """
        Remove unwanted files, e.g. non galaxies fro sample

        PARAMETER
            files_df
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

            results = pool.map(self._remove_fits_file, files_indexes)

        number_fail = sum(results)
        print(f"Fail with {number_fail} files")
    ###########################################################################
    def _remove_fits_file(self, file_index: "int")->"None":
        """
        Remove file
        PARAMETERS
            file_index: specobjid of the object
        """
        pass


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
            file_index: index of a galaxy in that galaxy that the
                frame passed to the constructor of the class

        OUTPUT
            0 for successful operation, 1 otherwise
        """

        file_row = files_df.loc[file_index]

        [file_directory, spectrum_name] = self._get_file_location(file_row)

        file_location = f"{file_directory}/{spectrum_name}.fits"

        if not super().file_exists(file_location, exit=False):
            print(file_location)

            return 1

        with counter.get_lock():
            counter.value += 1
            print(f"[{counter.value}] Get {spectrum_name}", end="\r")

        result = self._open_fits_file(
            file_index, file_location, spectrum_name
        )

        return result

    ###########################################################################
    def _open_fits_file(
        self,
        file_index: "int",
        file_location: "str",
        spectrum_name: "str",
        # number_attempts: "int" = 10
    ) -> "int":

        """
        Get and save wav, flux and ivar
        PARAMETER

            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor
            file_location: location of the fits file including extension
            spectrum_name: name of spectrum
                spec-{plate}-{mjd}-{fiberid}
            number_attempts: [no implement yet]
                how many times to try if warning

        OUTPUT
            0 for successful operation, 1 otherwise
        """

        save_to = f"{self.rest_frame_directory}/{spectrum_name}.npy"

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

    ###########################################################################
    def _get_file_location(self, file_row: "pd.row") -> "list":
        """
        PARAMETERS
            file_row:
        OUTPUTS
            return [file_directory, spectrum_name]
                file_directory: location of the spectrum fits file
                spectrum_name: f'spec-{plate}-{mjd}-{fiberid}'
        """


        [plate, mjd, fiberid, run2d] = self._galaxy_identifiers(file_row)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"{self.data_directory}/sas/dr16/sdss/spectro/redux"
            f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectrum_name]

    ###########################################################################
    def _galaxy_identifiers(self, file_row: "df.row") -> "list":
        """
        PARAMETER
            file_row : pd.row from the object data frame passed to
                the constructor of the class

        OUTPUT
            return [plate, mjd, fiberid, run2d]
                plate: self explanatory
                mjd: date
                fiberid: self explanatory
                run2d: PENDING

        """

        plate = f"{file_row['plate']:04}"
        mjd = f"{file_row['mjd']}"
        fiberid = f"{file_row['fiberid']:04}"
        run2d = f"{file_row['run2d']}"

        return [plate, mjd, fiberid, run2d]

###############################################################################
