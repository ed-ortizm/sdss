import ctypes
from functools import partial
import multiprocessing as mp
import os

####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

####################################################################
# from src.raw.worker import worker, init_worker

################################################################################
class RawData:
    """Handled raw data from sdss"""

    def __init__(
        self,
        galaxies_df: "pd.df",
        data_directory: "str",
        output_directory: "str",
        number_processes: "int",
    ):
        """
        PARAMETERS

            galaxies_df : data frame containing meta data from sdss
            data_directory : sdss raw data's directory
            output_directory : redshift corrected spectra's
                and meta data directory
            number_processes : number of processes to use with mp.Pool

        OUTPUT
            RawDataProcessing object
        """

        self.df = galaxies_df
        self.number_processes = number_processes

        self._check_directory(data_directory, exit=True)
        self.data_directory = data_directory

        self._check_directory(output_directory)
        self.data_output_directory = output_directory

        self._check_directory(f"{output_directory}/rest_frame")
        self.rest_frame_directory = f"{output_directory}/rest_frame"

    ###########################################################################
    def save_data_frame(self, name):

        location = f"{self.data_output_directory}/{name}"

        self.df.to_csv(path_or_buf=location, index=False)

    ###########################################################################
    def _add_columns_data_frame(self, data: "dictionary") -> "None":
        """
        Add class and subclass classification for galaxy data frame

        PARAMETERS
            data: dictionary with data to update data frame in place.
                data dictionary has the format
                data = {

                }

        """

        def get_lambda(idx):
            return lambda x: x[idx]

        f = get_lambda(0)
        indexes = list(map(f, data))

        f = get_lambda(1)
        classification = list(map(f, data))

        f = get_lambda(2)
        sub_class = list(map(f, data))

        self.df.loc[indexes, "class"] = classification
        self.df.loc[indexes, "subclass"] = sub_class

    ###########################################################################
    def get_raw_spectra(self):
        """
        Save data frame with all meta data
        """

        print(f"Saving raw redshift corrected spectra and meta-data!")

        # counter = mp.Value("i", 0)
        galaxy_indexes = range(self.df.shape[0])

        with mp.Pool(
            processes=self.number_processes,
            # initializer=init_worker,
            # initargs=(counter,),
        ) as pool:

            results = pool.map(self._get_spectra, galaxy_indexes)

        results = self._filter(results, parameter=1)

        self._add_columns_data_frame(results)

    ###########################################################################
    def _filter(self, results: "list", parameter) -> "list":
        """
        Filters input list

        PARAMETERS
            results: list from parallel computation.
                Visible satellites are a list, non visible is None

        OUTPUTS
            returns list with visible satellites
        """

        results = list(filter(lambda x: x != parameter, results))

        return results

    ###########################################################################
    def _get_spectra(self, galaxy_index: "int"):
        """
        Gets a row of meta data

        PARAMETERS
            galaxy_index: index of a galaxy in that galaxy that the
                frame passed to the constructor of the class

        OUTPUT
            meta data list, a row in the meta_data_frame:
                [
                name: galaxy name,
                z: redshift,
                snr: signal to noise ratio,
                run2d,
                sub-class: sdss classification,
                class: sdss classification
                ]
        """

        [file_directory, spectra_name] = self._get_file_location(galaxy_index)

        print(f"Process {spectra_name} --> N:{galaxy_index}", end="\r")

        file_location = f"{file_directory}/{spectra_name}.fits"

        if self._check_file(file_location, exit=False):

            return 1

        [classification, sub_class] = self._rest_frame(
            galaxy_index, file_location, spectra_name
        )

        meta_data = [galaxy_index, classification, sub_class]

        return meta_data

    ###########################################################################
    def _rest_frame(
        self, galaxy_index: "int", file_location: "str", spectra_name
    ):
        """
        De-redshifting

        PARAMETER
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor
            galaxy_file_location: directory location of the fits file

         OUTPUT
            return wave, flux, z, signal_noise_ratio, classification, sub_class
                wave: de-redshifted wave length
                flux: flux
                z: redshift
                signal_noise_ratio: signal to noise ratio
                classification: sdss pipeline classification
                sub_class: sdss subclass pipeline classification
        """

        with pyfits.open(file_location) as hdul:

            wave = 10.0 ** (hdul[1].data["loglam"])
            flux = hdul[1].data["flux"]
            classification = hdul[2].data["CLASS"][0]
            sub_class = hdul[2].data["SUBCLASS"][0]

        # deredshift flux?
        z = self.df.iloc[galaxy_index]["z"]
        z_factor = 1.0 / (1.0 + z)
        wave *= z_factor

        save_to = f"{self.rest_frame_directory}/{spectra_name}.npy"
        array_to_save = np.vstack((wave, flux))

        np.save(save_to, array_to_save)

        return [classification, sub_class]

    ###########################################################################
    def _get_file_location(self, galaxy_index: "int"):
        """
        PARAMETERS
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor

        OUTPUTS
            return sdss_directory, spectra_name, run2d

                sdss_directory: directory location of the spectra fits file
                spectra_name: f'spec-{plate}-{mjd}-{fiberid}'
                run2d: PENDING
        """

        galaxy = self.df.iloc[galaxy_index]
        [plate, mjd, fiberid, run2d] = self._galaxy_identifiers(galaxy)

        spectra_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"{self.data_directory}/sas/dr16/sdss/spectro/redux"
            f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectra_name]

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
    def _check_file(self, location: "str", exit: "bool" = False):
        """
        Check if a location exists, if not it creates it or
        exits depending on the value of exit
        """

        file_exists = not os.path.exists(location)

        if not file_exists and exit:

            print(f"File {location} NOT FOUND!")
            print("Code cannot execute")
            sys.exit()

        return file_exists


###############################################################################
