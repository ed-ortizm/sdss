#! /usr/bin/env python3

from argparse import ArgumentParser
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
ti = time.time()
################################################################################
parser = ArgumentParser()

parser.add_argument('--server', '-s', type=str)
parser.add_argument('--n_obs', type=int)
parser.add_argument('--normalization_type', '-n_type', type=str)

script_arguments = parser.parse_args()
################################################################################
n_obs = script_arguments.n_obs
normalization_type = script_arguments.normalization_type
local = script_arguments.server == 'local'
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

# Create and save interpolated data in master wave
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

# print(f'Sorting according to snMedian\n')
# SN_sorted_spectra = data_processing.sort_spec_SN(spectra=spectra)
################################################################################

print(f'Handling indefinite values')

spectra , wave = data_processing.indefinite_values_handler(spectra=spectra)

spectra = data_processing.missing_flux_replacement(spectra=spectra,
    method='median')

n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
print(f'Indefinite vals in the final array: {np.sum(n_indef)}')
###############################################################################
print(f'Normalizing data')
normalization_methods = ['median'] #, 'Z', 'min_max']

for method in normalization_methods:

    spectra = data_processing.normalize_spectra(spectra=spectra, method=method)
    ############################################################################
    print(f'Saving data')

    np.save(f'{spectra_path}/processed_spectra/spectra_{n_obs}_{method}.npy',
        spectra)
################################################################################
np.save(f'{spectra_path}/processed_spectra/wave_master_{n_obs}_processed.npy',
    wave)
################################################################################

tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')
