#! /usr/bin/env python3

import os
from time import time
import urllib
from glob import glob
from time import time, sleep
import multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

## From Guy Goren, Dovi's student

def getFitsFiles(gs, dbPath):

    """
    Downloads .fits files from SDSS dr16 for the galaxies specified in the gs DataFrame to dbPath (parallel version)
    :param gs: Pandas Dataframe: Containing the fields 'plate', 'mjd', 'fiberid', run2d'
    :param dbPath: String: Path to dB folder (where 'sas' folder is located).
    :param fitsFolder: String: fitsFolder name inside dbPath to place the .fits file
    :return: None.
    """

    print ('*** Getting  %d fits files ****' % (len(gs)))
    start = time()

    # Create fits' DB folder, if needed
    if not os.path.exists(dbPath):
        os.makedirs(dbPath)

    # Download .fits files parallelly
    pool = mp.Pool()
    f = partial(getFitFile_i, gs, dbPath)
    print('Starting first pool')
    res = pool.map(f, range(len(gs)))
    nFailed=sum(res)

    end = time()
    print('  Done! Finished downloading .fits files\t')
    print(f'  Failed to download: {nFailed} of {len(gs.index)}')
    print('  Time: ' + str(round(end - start, 2)) + ' [s]')


def getFitFile_i(gs, dbPath, i):
    """
    Downloading the ith .fits file in the gs DataFrame to dbPath.
    :param gs: Pandas Dataframe: Containing the fields 'plate', 'mjd', 'fiberid', run2d'
    :param dbPath: String: Path to dB folder (where 'sas' folder is located).
    :param i: Int: The index of the galaxy's .fits file to download
    :return: 0, if succeeded to download. 1 if failed.
    """

    # Extract identifiers
    obj = gs.iloc[i]
    plate = str(obj['plate']).zfill(4)
    mjd = str(obj['mjd'])
    fiberid = str(obj['fiberid']).zfill(4)
    run2d = str(obj['run2d'])

    # Try & Except a failed Download
    try:
        getFitsFile(plate, mjd, fiberid, run2d, dbPath)
        return 0
    except Exception as e:
        print ('Failed to DL: ' + 'https://data.sdss.org/' + 'sas/dr16/sdss/spectro/redux/' + run2d + '/spectra/lite/' + plate + '/' + '-'.join(['spec',plate,mjd,fiberid]) + '.fits')
        print(str(e))
        return 1


def getFitsFile(plate, mjd, fiberid, run2d, dbPath):

    """
    Download a .fits files by plate,mjd,fiberid & run2d to dBPath
    :param plate: String (4 digits): Plate number
    :param mjd: String: mjd
    :param fiberid: String (4 digits): fiber ID
    :param run2d: String: run2d
    :param dbPath: String: Path to dB folder (where 'sas' folder is located).
    :return: None.
    """

    # Calculate paths
    filename = '-'.join(['spec',plate,mjd,fiberid]) + '.fits'
    SDSSpath = '/sas/dr16/sdss/spectro/redux/' + run2d + '/spectra/lite/' + plate + '/'
    url = 'https://data.sdss.org/' + SDSSpath + filename
    folderPath = dbPath + SDSSpath
    dest = folderPath + filename

    # Create folder for the .fits file, if needed
    if not os.path.exists(folderPath):
        os.makedirs(folderPath)

    # Try to download the .fits file 10 times
    # A file under 60kB is considered a failure
    if not(os.path.isfile(dest)):
        urllib.request.urlretrieve(url, dest)

    j = 0
    while j < 10 and (os.path.getsize(dest) < 60000):
        os.remove(dest)
        urllib.request.urlretrieve(url, dest)
        j += 1
        sleep(1)
    if (os.path.getsize(dest) < 60000):
        print('Size is: ' + str(os.path.getsize(dest)) + ' '),
        os.remove(dest)
        raise Exception('Spectra wasn\'t found')

ti = time()

#### Obtain the spectra from SDSS DR16

## The sample
dbPath= f'/home/edgar/zorro/SDSSdata/'
gs = pd.read_csv(f'{dbPath}/gals_DR16.csv', header=0)

## Replacing z by z_noqso when possible
# z_noqso --> "Best redshift when excluding QSO fit in BOSS spectra
# (right redshift to use for galaxy targets)"

n0_rows = len(gs.index)
idx = (gs['z_noqso'] != 0)
n_z_noqso = len(gs.loc[idx,'z_noqso'].index)
#n_z_noqso = len(gs.loc[idx,'z_noqso'].shape[0])
print(f'There are {n0_rows} objects')
print(f'{n_z_noqso} redshifts were replaced by the z_noqso value')
print(f'{n0_rows - n_z_noqso} redshifts remained the same')

#gs.loc[idx,'z'] = gs.loc[idx,'z_noqso'].values
# This one is slighly faster
gs.loc[idx,'z'] = gs.z_noqso.loc[idx].values

# Remove galaxies with redshift z<=0.01

gs = gs[gs.z > 0.01]
n1_rows = len(gs.index)

print(f'{n0_rows - n1_rows} galaxies with z <= 0.01 removed')

gs.index = np.arange(n1_rows)

# Choose the top n_obs median SNR objects
gs.sort_values(by=['snMedian'], ascending=False, inplace=True)
gs.to_csv(f'{dbPath}/gs_SN_median_sorted.csv')

n_obs = 1_000 # 3188712
if not os.path.exists(f'{dbPath}/gs_{n_obs}.csv'):
    print(f'Creating file: gs_{n_obs}.csv')
    gs = gs[:n_obs]
    gs.index = np.arange(n_obs)

    # Create links for the Download

    url_head = 'http://skyserver.sdss.org/dr16/en/tools/explore/summary.aspx?plate='

    # It cannot be done with a big f string
    gs['url'] = url_head + gs['plate'].map('{:04}'.format) + '&mjd='\
                + gs['mjd'].astype(str) \
                + '&fiber=' + gs['fiberid'].map('{:04}'.format)
    gs.to_csv(f'{dbPath}/gs_{n_obs}.csv')
    tf = time()
    print(f'Building links for download took: {tf-ti:.2f} [s]')

    # Downloading the data...

    ti = time()

    getFitsFiles(gs,dbPath)

    tf = time()

    print(f'Data download took: {tf-ti:.2f} [s]')

# Plotting the z and SNR distribution

# fig, axarr = plt.subplots(1, 2)
#
# fig.set_figheight(7)
# fig.set_figwidth(14)
#
# ax1 = axarr[0]
# ax1.hist(gs.z, bins=int(len(gs)/50))
# ax1.set_xlabel('z')
# ax1.set_ylabel('Count')
# ax1.set_title('Redshift Histogram')
#
# ax2 = axarr[1]
# ax2.hist(gs.snMedian, bins=int(len(gs)/50))
# ax2.set_xlabel('median SNR')
# ax2.set_ylabel('Count')
# ax2.set_title('Median SNR Histogram')
#
# plt.show()
# plt.close
