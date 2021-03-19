#! /usr/bin/env python3
import glob
import os
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
# Loading data DataFrame with galaxies info
# SNR sorted data
gs = pd.read_csv(f'{spectra_path}/gals_DR16.csv')

# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)

# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(len(gs))

# Choose the top n_obs median SNR objects
n_obs = -1
if n_obs != -1:
    gs = gs[:n_obs]

################################################################################
# Data processing
data_processing = DataProcessing(galaxies_df= gs, n_processes=60)
data_processing.get_fluxes_SN()
################################################################################
# Getting array for train opening
dat_processing.
################################################################################

tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')
