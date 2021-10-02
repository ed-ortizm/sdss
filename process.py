#! /usr/bin/env python3
####################################################################
# Idea is to handle everything from the data frame
####################################################################
from configparser import ConfigParser, ExtendedInterpolation
import os
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
parser.read("process.ini")
################################################################################
# # Data processing
# Loading data DataFrame with galaxies info
data_directory = parser.get("directories", "data")
meta_data_file = parser.get("files", "meta_data")

galaxies_frame = pd.read_csv(f"{data_directory}/{meta_data_file}")
################################################################
number_processes = parser.getint("parameters", "processes")

data_process = data.DataProcess(
    galaxies_frame=galaxies_frame, number_processes=number_processes
)
################################################################
# interpolate spectra
wave_grid = data.get_grid(parser)

number_spectra = parser.getint("parameters", "spectra")

if number_spectra == -1:
    number_spectra = galaxies_frame.shape[0]

output_directory = parser.get("directories", "output")
output_directory = f"{output_directory}_{number_spectra}"

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

if not os.path.exists(f"{output_directory}/fluxes_interp.npy"):

    spectra = np.load(f"{output_directory}/fluxes_interp.npy")

else:

    spectra = data_process.interpolate(
        wave_master=wave_grid,
        data_directory=data_directory,
        output_directory=output_directory,
        number_spectra=number_spectra,
    )
################################################################
print(f"Handling indefinite values")

drop = parser.getfloat("parameters", "drop")

spectra, wave = data_process.drop_indefinite_values(
    spectra=spectra, wave_master=wave_grid, drop=drop
)

spectra = data_process.missing_flux_replacement(
    spectra=spectra, method="median"
)
################################################################################
print(f"Normalizing data")
#
spectra = data_process.normalize_spectra(spectra=spectra)
#
print(f"Saving data")
#
np.save(f"{output_directory}/fluxes.npy", spectra)
np.save(f"{output_directory}/wave.npy", wave)
################################################################################
t1 = time.time()
print(f"Running time: {t1-t0:.2f} [s]")
