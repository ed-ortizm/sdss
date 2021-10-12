import ctypes
from functools import partial
import multiprocessing as mp
import os

####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

################################################################################


class RawData:
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
    def _update_columns_data_frame(self, data_frame: "pandas.DataFrame"):
        """Add class and subclass classification for galaxy data frame"""
        # can I use a dictionary ?
        # d1 = {idx:clas, ...}, d2 = {idx:subclas}
        # df["class"] = d1
        # df["subclass"] = d2
        # Update in place
        pass

    ###########################################################################
    def get_raw_spectra(self):
        """
        Save data frame with all meta data
        """

        print(f"Saving raw redshift corrected spectra and meta-data!")

        galaxy_indexes = range(self.df.shape[0])

        with mp.Pool(processes=self.number_processes) as pool:
            results = pool.map(self._get_spectra, galaxy_indexes)

        self.meta_data_frame = pd.DataFrame(
            results,
            columns=[
                "galaxy_index",
                "name",
                "z",
                "snr",
                "run2d",
                "sub-class",
                "class",
            ],
        )

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

        [sdss_directory, spectra_name, run2d] = self._galaxy_localization(
            galaxy_index
        )

        print(f"Process {spectra_name} --> N:{galaxy_index}", end="\r")

        [plate, mjd, fiberid] = spectra_name.split("-")[1:]

        galaxy_fits_location = f"{sdss_directory}/{spectra_name}.fits"

        if not os.path.exists(galaxy_fits_location):

            print(f"[{galaxy_index}] NOT FOUND: {galaxy_fits_location}")

            meta_data = [
                galaxy_index,  # INDEX IN THE CURATED DATAFRAME
                spectra_name,
                np.nan,
                np.nan,
                run2d,
                np.nan,
                np.nan,
            ]

            return meta_data

        else:

            [
                wave,
                flux,
                z,
                signal_noise_ratio,
                classification,
                sub_class,
            ] = self._rest_frame(galaxy_index, galaxy_fits_location)

            np.save(
                f"{self.rest_frame_directory}/{spectra_name}.npy",
                np.vstack((wave, flux)),
            )

            meta_data = [
                galaxy_index,
                spectra_name,
                z,
                signal_noise_ratio,
                run2d,
                classification,
                sub_class,
            ]

            return meta_data

    ###########################################################################
    def _rest_frame(self, galaxy_index: "int", galaxy_fits_location: "str"):
        """
        De-redshifting

        PARAMETER
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor
            galaxy_fits_location: directory location of the fits file

         OUTPUT
            return wave, flux, z, signal_noise_ratio, classification, sub_class
                wave: de-redshifted wave length
                flux: flux
                z: redshift
                signal_noise_ratio: signal to noise ratio
                classification: sdss pipeline classification
                sub_class: sdss subclass pipeline classification
        """

        try:
            with pyfits.open(galaxy_fits_location) as hdul:

                wave = 10.0 ** (hdul[1].data["loglam"])
                flux = hdul[1].data["flux"]
                classification = hdul[2].data["CLASS"][0]
                sub_class = hdul[2].data["SUBCLASS"][0]

            z = self.df.iloc[galaxy_index]["z"]
            z_factor = 1.0 / (1.0 + z)
            wave *= z_factor

            signal_noise_ratio = self.df.iloc[galaxy_index]["snMedian"]

            return [
                wave,
                flux,
                z,
                signal_noise_ratio,
                classification,
                sub_class,
            ]

        except:
            wave = np.nan * np.empty(1)
            flux = np.nan * np.empty(1)
            z = np.nan
            signal_noise_ratio = np.nan
            classification = ""
            sub_class = ""

            return [
                wave,
                flux,
                z,
                signal_noise_ratio,
                classification,
                sub_class,
            ]

    ###########################################################################
    def _galaxy_localization(self, galaxy_index: "int"):
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
        plate, mjd, fiberid, run2d = self._galaxy_identifiers(galaxy)

        spectra_name = f"spec-{plate}-{mjd}-{fiberid}"

        sdss_directory = (
            f"{self.data_directory}/sas/dr16/sdss/spectro/redux"
            f"/{run2d}/spectra/lite/{plate}"
        )

        return [sdss_directory, spectra_name, run2d]

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

        return plate, mjd, fiberid, run2d

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


###############################################################################
