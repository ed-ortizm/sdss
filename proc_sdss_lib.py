import os
import urllib
from glob import glob
from time import time, sleep

import numpy as np
import astropy.io.fits as pyfits
import multiprocessing as mp
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from functools import partial
import matplotlib
import matplotlib.pyplot as plt

from constants_sdss import wave_master
from constants_sdss import working_dir, science_arxive_server_path
from constants_sdss import processed_spectra_path, spectra_path
## Me
################################################################################
def proc_spec(fnames):

    print('Processing all spectra')

    N = len(fnames)
    spec = np.empty((N, wave_master.size))

    for idx, fname in enumerate(fnames[:N]):
        print(f'Processing spectra N° {idx+1} --> {fname}', end='\r')
        spec[idx, :] = np.load(fname)

    print(f'indf vals: {np.count_nonzero(~np.isfinite(spec))}')

# Discarding spectrum with more than 10% of indefininte
# valunes in a given wl for al training set
    wkeep = np.where(np.count_nonzero(~np.isfinite(spec), axis=0) < spec.shape[0] / 10)
# Removing one dimensional axis since wkeep is a tuple
    spec = np.squeeze(spec[:, wkeep])
    wave_master = np.squeeze(spec[:, wkeep])

    print(f'indf vals: {np.count_nonzero(~np.isfinite(spec))}')

# Replacing indefinite values in a spectrum with its nan median
    for flx in spec.T:
        flx[np.where(~np.isfinite(flx))] = np.nanmedian(flx)

    print(f'indf vals: {np.count_nonzero(~np.isfinite(spec))}')

# Nomalize by the median and reduce noise with the standar deviation
    spec *= 1/np.median(spec, axis=1).reshape((spec.shape[0], 1))
#    spec *= 1/np.std(spec, axis=1).reshape((spec.shape[0], 1))


    np.save(f'spec_{N}.npy', spec)
    np.save(f'wave_master.npy', wave_master)

def get_spectra(gs, dbPath):
    """
    Computes the spectra interpolating over a master grid of wavelengths
    Parameters
    ----------
    gs : Pandas DataFrame with info of the galaxies
    dbPath : String : Path to data base

    Returns
    -------
    wave_master_grid : numpy array 1-D : The master wavelength grid
    flxs :  numpy array 2-D : Interpolated spectra over the grid
    """

    print(f'Getting grid of wavelengths and spectra from {len(gs)} .fits files')
    # http://python.omics.wiki/multiprocessing_map/multiprocessing_partial_function_multiple_arguments
    f = partial(flx_rest_frame_i, gs, dbPath)

    # close the pool (with) & do partial before
    with mp.Pool() as pool:
        pool.map(f, range(len(gs)))


#
    print('Job finished')


def flx_rest_frame_i(gs, dbPath, i):
    """
    Computes the min and max value in the wavelenght grid for the ith spectrum.
    Computes the interpolation functtion for the ith spectrum.
    Parameters
    ----------
    gs : Pandas DataFrame with info of the galaxies
    dbPath : String : Path to data base
    i : int : Index of .fits file in the DataFrame
    Returns
    -------
    min : float : minimum value if the wavelength grid for the ith spectrum
    max : float : maximun value if the wavelength grid for the ith spectrum
    flx_intp :  interp1d object : Interpolation function for the ith spectrum
    """

    obj = gs.iloc[i]
    plate = obj['plate']
    mjd = obj['mjd']
    fiberid = obj['fiberid']
    run2d = obj['run2d']
    z = obj['z']

    print(f'Processing spectrun N° {i}', end='\r')
    flx_rest_frame(plate, mjd, fiberid, run2d, z, dbPath)


def flx_rest_frame(plate, mjd, fiberid, run2d, z, dbPath):
    """
    Computes the min and max value in the wavelenght grid for the spectrum.
    Computes the interpolation functtion for the spectrum.
    Parameters
    ----------
    plate : int : plate number
    mjd : int : mjd of observation (days)
    fiberid : int : Fiber ID
    run2d : str : 2D Reduction version of spectrum
    z : float : redshift, replaced by z_noqso when available.
    (z_nqso --> Best redshift when excluding QSO fit in BOSS spectra (right redshift to use for galaxy targets))
    dbPath : String : Path to data base
    Returns
    -------
    wl_min : float : minimum value of the wavelength grid for the spectrum
    wl_max : float : maximun value of the wavelength grid for the spectrum
    flx_intp :  interp1d object : Interpolation function for the spectrum

    """
    # Path to the .fits file of the target spectrum
    fname = f'spec-{plate:04}-{mjd}-{fiberid:04}.fits'
    SDSSpath = f'{science_arxive_server_path}\
        dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate:04}'
    dest = f'{SDSSpath}/{fname}'


    if not(os.path.exists(dest)):
        print(f'File {fname} not found!')
        return None

    with pyfits.open(dest) as hdul:
        wl_rg = 10. ** (hdul[1].data['loglam'])
        flx = hdul[1].data['flux']


    # Deredshifting & min & max
    z_factor = 1./(1. + z)
    wl_rg *= z_factor
    flx = np.interp(wave_master, wl_rg, flx, left=np.nan, right=np.nan)

    np.save(
        f'{spectra_path}/{fname.split(".")[0]}_wave_master.npy',
        flx)
