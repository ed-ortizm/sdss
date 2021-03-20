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
################################################################################
ti = time.time()
################################################################################
# Loading data DataFrame with galaxies info
# SNR sorted data
gs = pd.read_csv(f'{spectra_path}/gals_DR16.csv')

# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)

# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(len(gs))

# Choose the top n_obs median SNR objects
if n_obs != -1:
    gs = gs[:n_obs]

################################################################################
# Data processing
data_processing = DataProcessing(galaxies_df= gs, n_processes=60)
data_processing.get_fluxes_SN()
################################################################################
# Getting array
fnames = glob.glob(
    f'{spectra_path}/interpolated_spectra/*interpolated.npy'
)
print(f'Number of files: {len(fnames)}')
spectra = data_processing.spec_to_single_array(fnames=fnames[:n_obs])
SN_sorted_spectra = data_processing.sort_spec_SN(spectra=spectra)
print(SN_sorted_spectra[:20, -1])
spectra , wave_master= data_processing.indefinite_values_handler(spectra=SN_sorted_spectra)
################################################################################

tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')
