#! /usr/bin/env python3
from glob import glob
import os
from time import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from proc_sdss_lib import get_spectra, proc_spec
################################################################################
working_directory = '/home/edgar/zorro/SDSSdata'
################################################################################
ti = time()
################################################################################
# Data processing

## Loading DataFrame with the data of the galaxies
# um I don't have SN_median sorted, got to get it
gs = pd.read_csv(f'{working_dir}/data/gs_SN_median_sorted.csv')


n_obs = 100_000 # 3188712
gs_n = gs[:n_obs]
gs_n.index = np.arange(n_obs)
get_spectra(gs_n, working_dir)

fnames = glob(f'{working_dir}/data/data_proc/*_wave_master.npy')

proc_spec(fnames[:])


tf = time()

print(f'Running time: {tf-ti:.2f} [seg]')
