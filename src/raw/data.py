import multiprocessing as mp
import os
import sys
import warnings

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

###############################################################################
def init_get_data_worker(input_counter: "mp.Value", in_df) -> "None":
    """
    Initialize worker for download
    PARAMETERS
        counter: counts the number of the child process
    """
    global counter
    global galaxies_df

    counter = input_counter
    galaxies_df = in_df

################################################################################
class GetRawData:
    """Handled raw data from sdss"""

    def __init__(
        self,
        data_directory: "str",
        output_directory: "str",
        number_processes: "int",
    ):
        """
        PARAMETERS

            galaxies_df : data frame containing meta data from sdss
            data_directory : sdss raw data's directory
            output_directory : save here .npy files
            number_processes : number of processes to use with mp.Pool

        OUTPUT
            GetRawData object
        """

        self.number_processes = number_processes

        self._check_directory(data_directory, exit=True)
        self.data_directory = data_directory

        self._check_directory(output_directory)
        self.data_output_directory = output_directory

        self._check_directory(f"{output_directory}/rest_frame")
        self.rest_frame_directory = f"{output_directory}/rest_frame"

    ###########################################################################
    def save_raw_data(self, galaxies_df: "pandas dataframe")->"None":
        """
        Save data frame with all meta dat

        PARAMETERS
            galaxies_df: data frame with meta data of galaxies
        """

        print(f"Save wave, flux, ivar!")

        galaxy_indexes = galaxies_df.index.values

        counter = mp.Value("i", 0)

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_get_data_worker,
            initargs=(counter,galaxies_df),
        ) as pool:

            results = pool.map(self._get_data, galaxy_indexes)


        number_fail = sum(results)
        print(f"Fail with {number_fail} files")

    ###########################################################################
    def _get_data(self, galaxy_index: "int")-> "int":
        """
        Gets a row of meta data

        PARAMETERS
            galaxy_index: index of a galaxy in that galaxy that the
                frame passed to the constructor of the class

        OUTPUT
            1 if file is not present, 0 otherwise
        """

        [
            file_directory,
            spectrum_name
        ] = self._get_file_location(galaxy_index)

        file_location = f"{file_directory}/{spectrum_name}.fits"

        if not self._file_exists(file_location, exit=False):
            print(file_location)

            return 1

        with counter.get_lock():
            counter.value += 1
            print(f"[{counter.value}] Get {spectrum_name}", end="\r")

        result = self._open_fits_file(
                                        galaxy_index,
                                        file_location,
                                        spectrum_name
                                    )

        return result

    ###########################################################################
    def _open_fits_file(
        self,
        galaxy_index: "int",
        file_location: "str",
        spectrum_name: "str",
        # number_attempts: "int" = 10
    )-> "None":
        """
        PARAMETER
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor
            file_location: directory location of the fits file
            spectrum_name:
            number_attempts:
        """

        save_to = f"{self.rest_frame_directory}/{spectrum_name}.npy"

        if self._file_exists(save_to, exit=False):
            print(f"Data of {spectrum_name} already saved!", end="\r")
            return 0

        warnings.filterwarnings(action="error")

        try:

            hdul = pyfits.open(file_location, memmap=False)

            wave = np.array(10.0 ** (hdul[1].data["loglam"]), dtype=float)
            flux = hdul[1].data["flux"]
            ivar = hdul[1].data['ivar']
            specobjid = int(hdul[2].data['specobjid'].item())

            df_specobjid = galaxies_df.iloc[galaxy_index]['specobjid']
            assert specobjid == df_specobjid , 'specobjid do not match'

            hdul.close()

            array_to_save = np.vstack((wave, flux, ivar))

            np.save(save_to, array_to_save)

            return 0

        except Exception as e:

            print(f"Problem with {spectrum_name}")
            print(e)

            return 1

    ###########################################################################
    def _get_file_location(self, galaxy_index: "int")-> "list":
        """
        PARAMETERS
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor

        OUTPUTS
            return [file_directory, spectrum_name]
                file_directory: location of the spectrum fits file
                spectrum_name: f'spec-{plate}-{mjd}-{fiberid}'
        """

        galaxy = galaxies_df.iloc[galaxy_index]

        [plate, mjd, fiberid, run2d] = self._galaxy_identifiers(galaxy)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"{self.data_directory}/sas/dr16/sdss/spectro/redux"
            f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectrum_name]

    ###########################################################################
    def _galaxy_identifiers(self, galaxy: "df.row"):
        """
        PARAMETER
            galaxy : pd.row from the galaxy data frame passed to
                the constructor of the class

        OUTPUT
            return plate, mjd, fiberid, run2d

                plate: self explanatory
                mjd: date
                fiberid: self explanatory
                run2d: PENDING

        """

        plate = f"{galaxy['plate']:04}"
        mjd = f"{galaxy['mjd']}"
        fiberid = f"{galaxy['fiberid']:04}"
        run2d = f"{galaxy['run2d']}"

        return [plate, mjd, fiberid, run2d]

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
    def _file_exists(self, location: "str", exit: "bool" = False):
        """
        Check if a location is a file, if not exits depending
        on the value of exit
        """

        file_exists = os.path.isfile(location)

        if not file_exists:

            file_name = location.split("/")[-1]

            if exit:
                print(f"File {file_name} NOT FOUND!")
                print("Code cannot execute")
                sys.exit()

            return file_exists

        return file_exists


###############################################################################
