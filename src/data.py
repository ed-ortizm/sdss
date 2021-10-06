from functools import partial
import multiprocessing as mp
import os

####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

################################################################################
class DataProcess:
    def __init__(
        self,
        galaxies_frame: "pd.df",
        number_processes: "int",
        grid_parameters: "dict",
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
    def interpolate(self, data_directory: "str", output_directory: "str"):
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory

        PARAMETERS
            data_directory:
            output_directory:

        OUTPUT
            mp.pool list with integers telling whether the process was successful or not.
            Interpolate spectra will be saved in output directory
        """
        print(f"Interpolate spectra...")

        number_spectra = self.frame.shape[0]
        fluxes = np.empty((number_spectra, self.grid.size))

        galaxy_names = self.frame.name
        spectra_directory = f"{data_directory}/rest_frame"
        self._check_directory(spectra_directory)

        for idx, galaxy_name in enumerate(galaxy_names):

            flux = self._interpolate(galaxy_name, spectra_directory)
            fluxes[idx, :] = flux[:]

        self._check_directory(output_directory)

        self.frame.to_csv(f"{output_directory}/meta_data.csv", index=False)

        np.save(f"{output_directory}/fluxes_interp.npy", fluxes)

        return fluxes

    ###########################################################################
    def _interpolate(self, galaxy_name, spectra_directory):

        print(f"Interpolate {galaxy_name}", end="\r")

        spectrum = np.load(f"{spectra_directory}/{galaxy_name}.npy")

        flux = np.interp(
            self.grid,
            spectrum[0],  # wave
            spectrum[1],  # flux
            left=np.nan,
            right=np.nan,
        )

        return flux
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


###############################################################################
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

        self.meta_data_frame = None

    ############################################################################
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

        with pyfits.open(galaxy_fits_location) as hdul:
            wave = 10.0 ** (hdul[1].data["loglam"])
            flux = hdul[1].data["flux"]
            classification = hdul[2].data["CLASS"][0]
            sub_class = hdul[2].data["SUBCLASS"][0]

        z = self.df.iloc[galaxy_index]["z"]
        z_factor = 1.0 / (1.0 + z)
        wave *= z_factor

        signal_noise_ratio = self.df.iloc[galaxy_index]["snMedian"]

        return wave, flux, z, signal_noise_ratio, classification, sub_class

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


###############################################################################
