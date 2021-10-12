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
from src.process import data

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

grid_parameters = dict(parser.items("grid"))

output_directory = parser.get("directories", "output")

data_process = data.DataProcess(
    galaxies_frame=galaxies_frame,
    number_processes=number_processes,
    grid_parameters=grid_parameters,
    data_directory=data_directory,
    output_directory=output_directory,
)
################################################################
# interpolate spectra
have_to_interpolate = parser.getboolean("parameters", "interpolate")

if have_to_interpolate:

    spectra = np.load(f"{output_directory}/fluxes_interp.npy")

else:
    
    spectra = data_process.interpolate()

################################################################
print(f"Handle indefinite values")

drop_fraction = parser.getfloat("parameters", "drop")

spectra, wave = data_process.drop_indefinite_values(
    spectra=spectra, drop=drop_fraction
)

print(f"Replace missing flux")

spectra = data_process.replace_missing_flux(spectra=spectra, method="median")
################################################################################
print(f"Normalize data")

spectra = data_process.normalize(spectra=spectra)

print(f"Save data")

np.save(f"{output_directory}/fluxes.npy", spectra)
np.save(f"{output_directory}/wave.npy", wave)
################################################################################
t1 = time.time()
print(f"Running time: {t1-t0:.2f} [s]")
