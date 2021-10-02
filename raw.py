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
# # Data processing
data_directory = parser.get("directories", "data")
gs_frame = parser.get("files", "raw_meta_data")

gs = pd.read_csv(f"{data_directory}/{gs_frame}")
####################################################################
# Use z_noqso when available
gs.z = np.where(gs.z_noqso.ne(0), gs.z_noqso, gs.z)
# Remove galaxies with redshift z<=0.01
gs = gs[gs.z > 0.01]
gs.index = np.arange(gs.shape[0])
####################################################################
number_spectra = parser.getint("parameters", "number_spectra")
if number_spectra != -1:
    gs = gs[:number_spectra]
####################################################################
output_directory = parser.get("directories", "output")
number_processes = parser.getint("parameters", "number_processes")

data_processing = data.RawData(
    galaxies_df=gs,
    data_directory=data_directory,
    output_directory=output_directory,
    number_processes=number_processes,
)
####################################################################
# get raw spectra
data_processing.get_raw_spectra()
################################################################################
# saving data frame with meta data of the raw a spectra in the rest frame
meta_data = parser.get("files", "meta_data")

data_processing.meta_data_frame.to_csv(
    path_or_buf=f"{output_directory}/{meta_data}", index=False
)
################################################################################
t1 = time.time()
print(f"Run time: {t1-t0}")
