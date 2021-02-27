#! /usr/bin/env python3
from glob import glob
import os
from time import time

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from constants_AEs import m_wl, data_proc
from proc_sdss_lib import get_spectra, proc_spec

ti = time()
# Data processing

## Loading DataFrame with the data of the galaxies

dbPath = f'/home/edgar/zorro/SDSSdata'
gs = pd.read_csv(f'{dbPath}/gs_SN_median_sorted.csv')


n_obs = 100_000 # 3188712
gs_n = gs[:n_obs]
gs_n.index = np.arange(n_obs)
get_spectra(gs_n, dbPath)

fnames = glob(f'{data_proc}/*.npy')

proc_spec(fnames[:])


tf = time()

print(f'Running time: {tf-ti:.2f} [seg]')
