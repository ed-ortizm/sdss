import os
import glob
import time
import urllib
import sys

import astropy.io.fits as pyfits
from functools import partial
import matplotlib
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import pandas as pd
from sklearn import preprocessing

from constants_sdss import n_waves, wave_master
from constants_sdss import working_dir, science_arxive_server_path
from constants_sdss import processed_spectra_path, spectra_path
################################################################################
class DataProcessing:

    def __init__(self, galaxies_df, n_processes): #fnames: list, SN_threshold: float):

        self.galaxies_df = galaxies_df
        self.n_processes = n_processes

        self.raw_spectra_path = f'{spectra_path}/raw_spectra'
        if not os.path.exists(self.raw_spectra_path):
            os.makedirs(self.raw_spectra_path, exist_ok=True)

        self.interpolated_spectra_path = f'{spectra_path}/interpolated_spectra'
        if not os.path.exists(self.interpolated_spectra_path):
            os.makedirs(self.interpolated_spectra_path, exist_ok=True)

        self.processed_spectra_path = f'{spectra_path}/processed_spectra'
        if not os.path.exists(self.processed_spectra_path):
            os.makedirs(self.processed_spectra_path, exist_ok=True)

    def normalize_spectra(self, spectra:'array', method:'str'):

        if method=='median':

            spectra[:, :-5] *= 1/np.median(spectra[:, :-5], axis=1).reshape(
                (spectra.shape[0], 1)
            )

        elif method=='Z':

            spectra[:, :-5] = preprocessing.scale(spectra[:, :-5])

        elif method=='min_max':

            spectra[:, :-5] = preprocessing.MinMaxScaler().fit_transform(
                spectra[:, :-5]
            )

        else:
            print(f'the only supported normalizations methods are:')
            print(f'median, Z (standarization) and min max normalization')
            sys.exit()

        return spectra

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

    def indefinite_values_handler(self, spectra: 'np.array',
        discard_fraction: 'float'=0.1):
        #global wave_master

        print(f'spectra shape before keep_spec_mask: {spectra.shape}')

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f'Indefinite vals in the input array: {np.sum(n_indef)}')

        keep_flux_mask =  n_indef < spectra.shape[0]*discard_fraction

        spectra = spectra[:, keep_flux_mask]
        print(f'spectra shape after keep_spec_mask: {spectra.shape}')

        wave = wave_master[keep_flux_mask[: -5]]

        n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
        print(f'Indefinite vals in the NEW array: {np.sum(n_indef)}')

        return spectra, wave

    def sort_spec_SN(self, spectra: 'array'):

        SN_arg_sort = np.argsort(-1*spectra[:, -1])
        spectra = spectra[SN_arg_sort, :]

        return spectra

    def spec_to_single_array(self, fnames: 'list'):

        n_spectra = len(fnames)
        n_fluxes = np.load(fnames[0]).size

        spectra = np.empty((n_spectra, n_fluxes))

        for idx, file_path in enumerate(fnames):

            fname = file_path.split('/')[-1].split('_')[0]

            print(f"Loading {fname} to single array", end='\r')

            spectra[idx, :] = np.load(file_path)

        return spectra

    def get_fluxes_SN(self):

        print(f'Saving raw and interpolated spectra')

        params = range(len(self.galaxies_df))

        with mp.Pool(processes=self.n_processes) as pool:
            res = pool.map(self._get_spectra, params)
            n_failed = sum(res)

        print(f'Spectra saved. Failed to save {n_failed}')

    def _get_spectra(self, idx_galaxy: int):

        galaxy_fits_path, fname = self._galaxy_fits_path(idx_galaxy)
        fname = fname.split('.')[0]
        [plate, mjd, fiberid] = fname.split('-')[1:]
        # print(f'plate: {plate}, mjd: {mjd}, fiberid: {fiberid}')

        if not os.path.exists(galaxy_fits_path):
            print(f'{fname} not found')
            return 1

        if os.path.exists(f'{self.interpolated_spectra_path}/{fname}.npy'):

            print(f'{fname} already processed', end='\r')
            return 0

        else:

            wave, flux, z, SN = self._rest_frame(idx_galaxy, galaxy_fits_path)
            flux_interpolated = np.interp(
                wave_master, wave, flux, left=np.nan, right=np.nan
            )

            np.save(f'{self.raw_spectra_path}/{fname}.npy',
                np.hstack(
                    (flux, int(plate), int(mjd), int(fiberid), z, SN)
                )
            )

            np.save(f'{self.interpolated_spectra_path}/{fname}_interpolated.npy',
                np.hstack(
                    (flux_interpolated, int(plate), int(mjd), int(fiberid), z, SN)
                )
            )

            return 0

    def _rest_frame(self, idx_galaxy, galaxy_fits_path):
        """De-redshifting"""

        with pyfits.open(galaxy_fits_path) as hdul:
            wave = 10. ** (hdul[1].data['loglam'])
            flux = hdul[1].data['flux']


        z = self.galaxies_df.iloc[idx_galaxy]['z']
        z_factor = 1./(1. + z)
        wave *= z_factor

        SN = self.galaxies_df.iloc[idx_galaxy]['snMedian']

        return wave, flux, z, SN

    def _galaxy_fits_path(self, idx_galaxy: int):

        galaxy = self.galaxies_df.iloc[idx_galaxy]
        plate, mjd, fiberid, run2d = self._galaxy_identifiers(galaxy)

        fname = f'spec-{plate}-{mjd}-{fiberid}.fits'

        SDSSpath = f'sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}'
        retrieve_path = f'{spectra_path}/{SDSSpath}'

        return f'{retrieve_path}/{fname}', fname

    def _galaxy_identifiers(self, galaxy):

        plate = f"{galaxy['plate']:04}"
        mjd = f"{galaxy['mjd']}"
        fiberid = f"{galaxy['fiberid']:04}"
        run2d = f"{galaxy['run2d']}"

        return plate, mjd, fiberid, run2d
################################################################################
################################################################################
class DownloadData:

    def __init__(self, files_data_frame, download_path, n_processes):
        """
        files_data_frame: Pandas DataFrame with all the imformation of the sdss
        galaxies

        download_path: (string) Path where the data will be downloaded
        """
        self.data_frame = files_data_frame
        self.download_path = download_path
        self.n_processes = n_processes

    def get_files(self):

        print(f'*** Getting {len(self.data_frame)} fits files ****')
        start_time_download = time.time()

        if not os.path.exists(self.download_path):
            os.makedirs(self.download_path)

        params = range(len(self.data_frame))

        with mp.Pool(processes=self.n_processes) as pool:
            res = pool.map(self._get_file, params)
            n_failed = sum(res)

        finish_time_download = time.time()

        print(f'Done! Finished downloading .fits files...')
        print(f'Failed to download {n_failed} files' )
        print(
            f'Download took {finish_time_download-start_time_download:.2f}[s]')


    def _get_file(self, idx_data_frame):

        object = self.data_frame.iloc[idx_data_frame]
        plate, mjd, fiberid, run2d = self._file_identifier(object)

        fname = f'spec-{plate}-{mjd}-{fiberid}.fits'

        SDSSpath = f'sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}'
        folder_path = f'{self.download_path}/{SDSSpath}'

        url =\
        f'https://data.sdss.org/{SDSSpath}/{fname}'

        if not os.path.exists(folder_path):
            os.makedirs(folder_path, exist_ok=True)

        # Try & Except a failed Download

        try:
            self._retrieve_url(url, folder_path, fname)
            return 0

        except Exception as e:

            print(f'Failed : {url}')

            print(f'{e}')
            return 1

    def _retrieve_url(self, url, folder_path, fname):

        if not(os.path.isfile(f'{folder_path}/{fname}')):

            print(f'Downloading {fname}', end='\r')

            urllib.request.urlretrieve(url, f'{folder_path}/{fname}')

            file_size = os.path.getsize(f'{folder_path}/{fname}')

            j = 0

            while j < 10 and (file_size < 60000):

                os.remove(f'{folder_path}/{fname}')
                urllib.request.urlretrieve(url, f'{folder_path}/{fname}')
                j += 1
                time.sleep(1)

            file_size = os.path.getsize(f'{folder_path}/{fname}')

            if file_size < 60000:
                print(f"Size of {fname}: {file_size}... Removing file!!")
                os.remove(f'{folder_path}/{fname}')
                raise Exception("Spectra wasn't found")

        else:
            print(f'{fname} already downloaded!!', end='\r')


    def _file_identifier(self, object):

        plate = f"{object['plate']:04}"
        mjd = f"{object['mjd']}"
        fiberid = f"{object['fiberid']:04}"
        run2d = f"{object['run2d']}"

        return plate, mjd, fiberid, run2d
