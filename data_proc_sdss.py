#! /usr/bin/env python3
from glob import glob
import os
from time import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from constants_sdss import science_arxive_server_path, spectra_path
from proc_sdss_lib import DownloadData, get_spectra, proc_spec
################################################################################
working_directory = '/home/edgar/zorro/spectra_sdss'
################################################################################
ti = time()
################################################################################
# Sort by SNR
gs = pd.read_csv(f'{spectra_path}/gals_DR16.csv')

# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)
#gs['z'] = [row['z_noqso'] if row['z_noqso']!=0 else row['z'] for i, row in gs.iterrows()]
#gs['test'] = [row['z_noqso'] if row['z_noqso']!=0 else row['z'] for i, row in gs.iterrows()]
#n = np.count_nonzero(np.where(gs.z==gs.test, True, False))
# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(len(gs))
#
# # Choose the top n_obs median SNR objects
n_obs = -1

if n_obs != -1:
    gs = gs[:n_obs]

# Data Download

download_spectra = DownloadData(
    files_data_frame=gs, download_path=spectra_path, n_processes=50)

# # Data processing
#
# ## Loading DataFrame with the data of the galaxies
# # um I don't have SN_median sorted, got to get it
# gs = pd.read_csv(f'{working_dir}/data/gs_SN_median_sorted.csv')
#
#
# gs_n.index = np.arange(n_obs)
# get_spectra(gs, spectra_path)
#
#fnames = glob(f'{working_dir}/data/data_proc/*_wave_.npy')
#
# proc_spec(fnames[:])


tf = time()

print(f'Running time: {tf-ti:.2f} [seg]')
