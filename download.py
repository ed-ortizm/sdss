#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import time
##############################################################
import numpy as np
import pandas as pd
##############################################################
from src import download
################################################################################
ti = time.time()
################################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read('download.ini')
################################################################################
# SNR sorted data
data_directory = parser.get('directories', 'data')
data_file_name = parser.get('files', 'data')
##############################################################
gs = pd.read_csv(f'{data_directory}/{data_file_name}')
# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)
# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(gs.shape[0])
##############################################################
number_spectra = parser.getint('parameters', 'number_spectra')
if number_spectra != -1:
    gs = gs[:number_spectra]
##############################################################
# Data Download
number_processes = parser.getint('parameters', 'number_processes')

download_spectra = download.DownloadData(
    files_data_frame=gs,
    download_path=data_directory,
    n_processes=number_processes
    )

download_spectra.get_files()
################################################################################
tf = time.time()

print(f'Running time: {tf-ti:.2f} [seg]')
