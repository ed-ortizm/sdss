#! /usr/bin/env python3
import glob
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from constants_sdss import spectra_path
from lib_processing_sdss import DataProcessing

################################################################################
n_obs = int(sys.argv[1])
local = sys.argv[2] == 'local'
################################################################################
ti = time.time()
################################################################################
# Loading data DataFrame with galaxies info
# SNR sorted data
if local:
    gs = pd.read_csv(f'{spectra_path}/gals_DR16_1000.csv')
else:
    gs = pd.read_csv(f'{spectra_path}/gals_DR16.csv')

# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)

# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(len(gs))

# Choose the top n_obs median SNR objects
if not local:
    if n_obs != -1:
        gs = gs[:n_obs]

################################################################################
# Data processing
data_processing = DataProcessing(galaxies_df= gs, n_processes=60)

if not local:
    data_processing.get_fluxes_SN()
################################################################################
# Getting array
fnames = glob.glob(
    f'{spectra_path}/interpolated_spectra/*_interpolated.npy'
)
if local:
    n_obs = len(fnames)

print(f'Number of files: {n_obs}')

print(f'spec to single array')
spectra = data_processing.spec_to_single_array(fnames=fnames[:n_obs])

print(f'Sorting according to snMedian')
SN_sorted_spectra = data_processing.sort_spec_SN(spectra=spectra)

print(f'Handling indefinite values')
spectra , wave = data_processing.indefinite_values_handler(
    spectra=SN_sorted_spectra
)

np.save(f'spectra_{n_obs}.npy', spectra)
np.save(f'wave_master_processed.npy', wave)
################################################################################

tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')
