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
parser.read("raw.ini")
################################################################################
# Data processing
raw_frame_location = parser.get("files", "raw_meta_data")

galaxy = pd.read_csv(raw_frame_location)
####################################################################
# Use z_noqso when available
galaxy.z = np.where(galaxy.z_noqso.ne(0), galaxy.z_noqso, galaxy.z)
# Remove galaxies with redshift z<=0.01
galaxy = galaxy[galaxy.z > 0.01]
galaxy.index = np.arange(galaxy.shape[0])
####################################################################
number_spectra = parser.getint("parameters", "number_spectra")
if number_spectra != -1:
    galaxy = galaxy[:number_spectra]
####################################################################
data_directory = parser.get("directories", "data")
output_directory = parser.get("directories", "output")
number_processes = parser.getint("parameters", "number_processes")

data_processing = data.RawData(
    galaxies_df=galaxy,
    data_directory=data_directory,
    output_directory=output_directory,
    number_processes=number_processes,
)
####################################################################
# get raw spectra
print("Get raw spectra")
data_processing.get_raw_spectra()
################################################################################
# saving data frame with meta data of the raw a spectra in the rest frame
meta_data_location = parser.get("files", "meta_data")

data_processing.meta_data_frame.to_csv(
    path_or_buf=meta_data_location,
    index=False
)
################################################################################
t1 = time.time()
print(f"Run time: {t1-t0}")
