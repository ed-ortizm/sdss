#! /usr/bin/env python3

#from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation

import numpy as np
import pandas as pd

from lib_processing_sdss import RawData
################################################################################
# parser = ArgumentParser()
# parser.add_argument('--server', '-s', type=str)
# script_arguments = parser.parse_args()
# local = script_arguments.server == 'local'
################################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read('config.ini')
# relevant directory
data_directory = parser.get('directories', 'data_directory')
output_directory = parser.get('directories', 'output_directory')

number_spectra = parser.getint('script parameters', 'number_spectra')
number_processes = parser.getint('script parameters', 'number_processes')
################################################################################
# Data processing
print(f'{data_directory}')
gs = pd.read_csv(f'{data_directory}/gals_DR16.csv')

# Use z_noqso if possible
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)

# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(len(gs))

if number_spectra != -1:
    gs = gs[:number_spectra]
################################################################################
data_processing = RawData(galaxies_df=gs,
    data_directory=data_directory, output_directory=output_directory,
    number_processes=number_processes)
################################################################################
#get raw spectra
data_processing.get_raw_spectra()
################################################################################
# saving data frame with meta data of the raw a spectra in the rest frame
print(data_processing.meta_data_frame)

data_processing.meta_data_frame.to_csv(
    path_or_buf=f'{output_directory}/meta_data_frame.csv')
