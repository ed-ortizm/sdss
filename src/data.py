from functools import partial
import multiprocessing as mp
import os
####################################################################
import astropy.io.fits as pyfits
import numpy as np
import pandas as pd
################################################################################
def get_grid(
    parser:'ConfigtParser obj'
    )->'np.array':
    """
    Computes the master grid for the interpolation of the spectra

    ARGUMENTS

        parser: ConfigurationParser object that contains information
        for this computation read from the .ini file

    RETURN
        wave_grid:
    """
    number_waves = parser.getint('constants', 'wave_master')
    master_lower = parser.getint('constants', 'master_lower')
    master_upper = parser.getint('constants', 'master_upper')

    wave_grid = np.linspace(
        master_lower,
        master_upper,
        number_waves
        )

    return wave_grid
################################################################################
class DataProcess:

    def __init__(self,
        galaxies_frame:'pd.df',
        number_processes:'int'
        ):
        """
        Class to process rest frame spectra
        INPUTS
            galaxy_frame: meta data of sdss galaxies, such as name,
                signal to noise ratio and z
            number_processes: number of jobs when processing a bulk of a spectra
        OUTPUT
            check how to document the constructor of a class
        """
        self.frame = galaxies_frame
        self.number_processes = number_processes
        self.fluxes = None

        # Single interpolated array,
        # count pythonic number of files in data diretory
    ############################################################################
    def interpolate(self,
        wave_master:'np.array',
        data_directory:'str',
        output_directory:'str',
        number_spectra:'int'=False
        ):
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory

        INPUTS
            wave_master: 1 dimensional array containing the common grid
                to use with all spectra
            data_directory:
            output_directory:
            number_spectra: number of spectra to process when testing

        OUTPUT
            mp.pool list with integers telling whether the process was successful or not.
            Interpolate spectra will be saved in output directory
        """
        print(f'Interpolate spectra...')

        # use a partial from itertools for interpolate function
        if number_spectra:
            index_galaxies = range(number_spectra)
        else:
            index_galaxies = range(self.frame.shape[0])

        self.fluxes = np.empty(len(index_galaxies), wave_master.size)

        interpolator = partial(
            self.interpolate_single,
            wave_master=wave_master,
            data_directory=data_directory,
            output_directory=output_directory
            )

        with mp.Pool(processes=self.number_processes) as pool:
            results = pool.map(interpolator, index_galaxies)
            number_failed = sum(results)

        print(f'Spectra saved. Failed to save {number_failed}')

        # np.save() save here in train dir self.fluxes and new metadataframe as well

    def interpolate_single(self,
        galaxy_index:'int',
        wave_master:'np.array',
        data_directory:'str',
        output_directory:'str'
        ):
        """
        Function to interpolate single spectrum to wav master

        INPUT
            galaxy_name: index of galaxy in the meta data frame
            wave_master: 1 dimensional array containing the common grid
                to use with all spectra
            data_directory:
            output_directory:
        OUTPUT
            interpolated spectrum as a numpy array
        """
        galaxy_name = self.frame.name[galaxy_index]
        spectrum_direction = f'{data_directory}/rest_frame/{galaxy_name}.npy'

        if os.path.exists(spectrum_direction):

            spectrum = np.load(spectrum_direction)

        else:

            print(f'There is no file: {galaxy_name}')

            return 1

        self.fluxes[galaxy_index, :] = np.interp(
            wave_master,
            spectrum[0], # wave
            spectrum[1], # flux
            left=np.nan,
            right=np.nan
            )
        # flux_direction = f'{output_directory}/{galaxy_name}.npy'
        # np.save(flux_direction, flux)
        return 0

    def normalize_spectra(self, spectra:'np.array'):
        # """"""

        # spectra[:, :] *= 1/np.median(spectra[:, :], axis=1).reshape(
        #         (spectra.shape[0], 1))

       #  return spectra
       pass

    def missing_flux_replacement(self, spectra:'array', method:'str'='median'):

        if method=='median':

            mask_replacement = ~np.isfinite(spectra)
            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmedian(spectra[idx, :])

        elif method=='mean':

            mask_replacement = ~np.isfinite(spectra)
            for idx, mask in enumerate(mask_replacement):
                spectra[idx, mask] = np.nanmean(spectra[idx, :])

        return spectra

    def indefinite_values_handler(self, spectra:'np.array',
        discard_fraction:'float'=0.1):
        #global wave_master

        print(f'spectra shape before keep_spec_mask: {spectra.shape}')

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f'Indefinite vals in the input array: {np.sum(n_indef)}')

        keep_flux_mask =  n_indef < spectra.shape[0]*discard_fraction

        spectra = spectra[:, keep_flux_mask]
        print(f'spectra shape after keep_spec_mask: {spectra.shape}')

        wave = wave_master[keep_flux_mask[: -8]]

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f'Indefinite vals in the NEW array: {np.sum(n_indef)}')

        return spectra, wave

    def spec_to_single_array(self, fnames: 'list'):

        n_spectra = len(fnames)
        n_fluxes = np.load(fnames[0]).size

        spectra = np.empty((n_spectra, n_fluxes))

        for idx, file_path in enumerate(fnames):

            fname = file_path.split('/')[-1].split('_')[0]

            print(f"Loading {fname} to single array", end='\r')

            spectra[idx, :] = np.load(file_path)

        return spectra
################################################################################
class RawData:

    def __init__(self, galaxies_df:'pd.df',
        data_directory:'str', output_directory:'str',
        number_processes:'int'):
        """
        INPUT

        galaxies_df : data frame containing meta data from sdss spectra
        data_directory : sdss raw data's directory
        output_directory : redshift corrected spectra's and meta data directory
        number_processes : number of processes to use with mp.Pool

        OUTPUT
        RawDataProcessing object
        """

        self.df = galaxies_df
        self.number_processes = number_processes

        self.data_directory = data_directory

        if not os.path.exists(self.data_directory):
            print(f'Path: {self.data_output_directory} does not exist!')

        self.data_output_directory = output_directory
        self.rest_frame_directory = f"{output_directory}/rest_frame"

        if not os.path.exists(self.data_output_directory):
            os.makedirs(self.data_output_directory)
            os.makedirs(f"{self.rest_frame_directory}")


        self.meta_data_frame = None
    ############################################################################
    def get_raw_spectra(self):
        """
        Save data frame with all meta data

        INPUTS
            None

        OUTPUT
            None
        """

        print(f'Saving raw redshift corrected spectra and meta-data!')

        galaxy_indexes = range(self.df.shape[0])

        with mp.Pool(processes=self.number_processes) as pool:
            results = pool.map(self._get_spectra, galaxy_indexes)

        self.meta_data_frame = pd.DataFrame(
            results,
            columns=['name', 'z', 'snr', 'run2d', 'sub-class', 'class'])

    def _get_spectra(self, galaxy_index:'int'):
        """
        Gets a row of meta data

        INPUT
            galaxy_index: index of a galaxy in that galaxy that the frame
            passed to the constructor of the class

        OUTPUT
            meta data list, intended to be a row in the meta_data_frame:
                meta_data_frame: [name: galaxy name,
                z: redshift,
                snr: signal to noise ratio,
                run2d,
                sub-class: sdss classification,
                class: sdss classification]
        """

        sdss_directory, spectra_name, run2d = self._galaxy_localization(
            galaxy_index)

        [plate, mjd, fiberid] = spectra_name.split('-')[1:]

        galaxy_fits_location = f'{sdss_directory}/{spectra_name}.fits'

        if not os.path.exists(galaxy_fits_location):

            print(f'{spectra_name}.fits file not found')

            meta_data = [spectra_name, np.nan, np.nan, run2d, np.nan, np.nan]

            return meta_data

        else:

            [wave, flux, z, signal_noise_ratio,
                classification, sub_class] = self._rest_frame(galaxy_index,
                                                 galaxy_fits_location)

            np.save(f"{self.rest_frame_directory}/{spectra_name}.npy",
                np.vstack((wave, flux)))

            meta_data = [spectra_name, z, signal_noise_ratio, run2d,
                classification, sub_class]

            return meta_data

    def _rest_frame(self, galaxy_index:'int', galaxy_fits_location:'str'):
        """
        De-redshifting

        INPUT
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
            wave = 10. ** (hdul[1].data['loglam'])
            flux = hdul[1].data['flux']
            classification = hdul[2].data['CLASS'][0]
            sub_class = hdul[2].data['SUBCLASS'][0]


        z = self.df.iloc[galaxy_index]['z']
        z_factor = 1./(1. + z)
        wave *= z_factor

        signal_noise_ratio = self.df.iloc[galaxy_index]['snMedian']

        return wave, flux, z, signal_noise_ratio, classification, sub_class

    def _galaxy_localization(self, galaxy_index:'int'):
        """
        INPUTS
            galaxy_index: index of the galaxy in the input data frame
                passed to the constructor

        OUTPUTS
            return sdss_directory, spectra_name, run2d

                sdss_directory: directory location of the spectra fits file
                spectra_name: f'spec-{plate}-{mjd}-{fiberid}'
                run2d: PENDING
        """

        galaxy = self.df.iloc[galaxy_index]
        plate, mjd, fiberid, run2d = self.galaxy_identifiers(galaxy)

        spectra_name = f'spec-{plate}-{mjd}-{fiberid}'

        sdss_directory = (f'{self.data_directory}/sas/dr16/sdss/spectro/redux'
            f'/{run2d}/spectra/lite/{plate}')

        return sdss_directory, spectra_name, run2d
    ############################################################################
    def galaxy_identifiers(self, galaxy:'df.row'):
        """
        INPUT
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
################################################################################
