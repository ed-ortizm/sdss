from functools import partial
import os
import sys
import time
import urllib

import astropy.io.fits as pyfits
import multiprocessing as mp
import numpy as np
import pandas as pd
################################################################################
class FitsPath:
    ############################################################################
    def __init__(self, galaxies_df, n_processes):

        self.galaxies_df = galaxies_df
        self.n_processes = n_processes
    ############################################################################
    def decode_base36(self, sub_class:'str'):

        if ' ' in sub_class:
            sub_class = sub_class.replace(' ', '')

        elif sub_class == '':
            sub_class = 'EC'

        return int(sub_class, 36)
    ############################################################################
    def encode_base36(self, sub_class:'int'):

        alphabet, base36 = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ', '']

        sub_class = abs(int(sub_class))

        while sub_class:
            sub_class, i = divmod(sub_class, 36)
            base36 = alphabet[i] + base36

        return base36 or alphabet[0]
    ############################################################################
    def class_sub_class(self, fits_paths:'list'):

        with mp.Pool(processes=self.n_processes) as pool:
            res = pool.map(self._class_sub_class, fits_paths)

        return res
    ############################################################################
    def _class_sub_class(self, fits_path:'str'):

        with pyfits.open(fits_path) as hdul:
            classification = hdul[2].data['CLASS']
            sub_class = hdul[2].data['SUBCLASS']

        return sub_class #[classification, sub_class]
    ############################################################################
    def get_all_paths(self):

        params = range(len(self.galaxies_df))

        with mp.Pool(processes=self.n_processes) as pool:
            res = pool.map(self._get_path, params)

        return res
    ############################################################################
    def _get_path(self, galaxy_index:'int'):

        galaxy_fits_path, fname = self._galaxy_fits_path(galaxy_index)

        return galaxy_fits_path
    ############################################################################
    def _galaxy_fits_path(self, galaxy_index:'int'):

        galaxy = self.galaxies_df.iloc[galaxy_index]
        plate, mjd, fiberid, run2d = self._galaxy_identifiers(galaxy)

        fname = f'spec-{plate}-{mjd}-{fiberid}.fits'

        SDSSpath = f'sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}'
        retrieve_path = f'{spectra_path}/{SDSSpath}'

        return f'{retrieve_path}/{fname}', fname
    ############################################################################
    def _galaxy_identifiers(self, galaxy):

        plate = f"{galaxy['plate']:04}"
        mjd = f"{galaxy['mjd']}"
        fiberid = f"{galaxy['fiberid']:04}"
        run2d = f"{galaxy['run2d']}"

        return plate, mjd, fiberid, run2d
################################################################################
