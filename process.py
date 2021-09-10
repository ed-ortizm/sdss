#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import time
####################################################################
import numpy as np
import pandas as pd
####################################################################
from src import data
################################################################################
t0 = time.time()
################################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read('process.ini')
################################################################################
# # Data processing
# Loading data DataFrame with galaxies info
data_directory = parser.get('directories', 'data')
meta_data_file = parser.get('files', 'meta_data')

galaxies_frame = pd.read_csv(f'{data_directory}/{meta_data_file}')
################################################################
number_processes = parser.getint('parameters', 'processes')

data_process = data.DataProcess(
    galaxies_frame=galaxies_frame,
    number_processes=number_processes)
################################################################
# interpolate spectra
wave_grid = data.get_grid(parser)
output_directory = parser.get('directories', 'output')
data_process.interpolate(
    wave_master=wave_grid,
    data_directory=data_directory,
    output_directory=output_directory
    )
################################################################
# spectra_name = galaxy_frame.name[:]
# # ################################################################################
# # # Getting array
# # fnames = glob.glob(
# #     f'{spectra_path}/interpolated_spectra/*_interpolated.npy'
# # )
# # print(f'Number of files: {n_obs}')
# #
# # print(f'spec to single array')
# # spectra = data_processing.spec_to_single_array(fnames=fnames[:n_obs])
# #
# # # print(f'Sorting according to snMedian\n')
# # # SN_sorted_spectra = data_processing.sort_spec_SN(spectra=spectra)
# # ################################################################################
# #
# # print(f'Handling indefinite values')
# #
# # spectra , wave = data_processing.indefinite_values_handler(spectra=spectra)
# #
# # spectra = data_processing.missing_flux_replacement(spectra=spectra,
# #     method='median')
# #
# # n_indef = np.count_nonzero(~np.isfinite(spectra), axis=0)
# # print(f'Indefinite vals in the final array: {np.sum(n_indef)}')
# # ###############################################################################
# # print(f'Normalizing data')
# # normalization_methods = ['median'] #, 'Z', 'min_max']
# #
# # for method in normalization_methods:
# #
# #     spectra = data_processing.normalize_spectra(spectra=spectra, method=method)
# #     ############################################################################
# #     print(f'Saving data')
# #
# #     np.save(f'{spectra_path}/processed_spectra/spectra_{n_obs}_{method}.npy',
# #         spectra)
# # ################################################################################
# # np.save(f'{spectra_path}/processed_spectra/wave_master_{n_obs}_processed.npy',
# #     wave)
# ################################################################################
t1 = time.time()
print(f'Running time: {t1-t0:.2f} [s]')
