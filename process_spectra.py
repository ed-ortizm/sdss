#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import time

import numpy as np
import pandas as pd

from lib_processing_sdss import DataProcess
################################################################################
t0 = time.time()
################################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read('process.ini')
    ############################################################################
# relevant directory and files
data_directory = parser.get('directories', 'data')
meta_data_file = parser.get('files', 'meta_data')
output_directory = parser.get('directories', 'output')
    ############################################################################
# script parameters
number_spectra = parser.getint('parameters', 'number_spectra')
number_processes = parser.getint('parameters', 'number_processes')

# wavelength master array
wave_master = parser.getint('constants', 'wave_master')
wave_master_lower = parser.getint('constants', 'wave_master_lower')
wave_master_upper = parser.getint('constants', 'wave_master_upper')
wave_master = np.linspace(wave_master_lower, wave_master_upper, wave_master)
################################################################################
# Loading data DataFrame with galaxies info
galaxies_frame = pd.read_csv(f'{data_directory}/{meta_data_file}')

# Data processing
data_process = DataProcess(
    galaxies_frame=galaxies_frame,
    number_processes=number_processes)

# interpolate spectra
# spectra_name = galaxy_frame.name[:]
data_process.interpolate(wave_master=wave_master,
                        data_directory=data_directory,
                        output_directory=output_directory)
# ################################################################################
# # Getting array
# fnames = glob.glob(
#     f'{spectra_path}/interpolated_spectra/*_interpolated.npy'
# )
# print(f'Number of files: {n_obs}')
#
# print(f'spec to single array')
# spectra = data_processing.spec_to_single_array(fnames=fnames[:n_obs])
#
# # print(f'Sorting according to snMedian\n')
# # SN_sorted_spectra = data_processing.sort_spec_SN(spectra=spectra)
# ################################################################################
#
# print(f'Handling indefinite values')
#
# spectra , wave = data_processing.indefinite_values_handler(spectra=spectra)
#
# spectra = data_processing.missing_flux_replacement(spectra=spectra,
#     method='median')
#
# n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
# print(f'Indefinite vals in the final array: {np.sum(n_indef)}')
# ###############################################################################
# print(f'Normalizing data')
# normalization_methods = ['median'] #, 'Z', 'min_max']
#
# for method in normalization_methods:
#
#     spectra = data_processing.normalize_spectra(spectra=spectra, method=method)
#     ############################################################################
#     print(f'Saving data')
#
#     np.save(f'{spectra_path}/processed_spectra/spectra_{n_obs}_{method}.npy',
#         spectra)
# ################################################################################
# np.save(f'{spectra_path}/processed_spectra/wave_master_{n_obs}_processed.npy',
#     wave)
################################################################################
t1 = time.time()
print(f'Running time: {t1-t0:.2f} [s]')
