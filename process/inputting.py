"""Imputting of missing values in interpolated spectra"""
from configparser import ConfigParser, ExtendedInterpolation
import time
from matplotlib.pyplot import grid

import numpy as np
import pandas as pd

from sdss.process import inputting
from sdss.utils.configfile import ConfigurationFile


start_time = time.time()

parser = ConfigParser(interpolation=ExtendedInterpolation())
name_config_file = "inputting.ini"
parser.read(f"{name_config_file}")

# A load data frame with meta data
meta_data_directory = parser.get("directory", "meta_data")

spectra_df_name = parser.get("files", "spectra_df")
spectra_df = pd.read_csv(
    f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
)

data_directory = parser.get("directory", "data")
# Load interpolated spectra
spectra_file_name = parser.get("files", "spectra")
spectra = np.load(f"{data_directory}/{spectra_file_name}")
# Load indexes and specobjid of interpolated spectra
ids_file_name = parser.get("files", "ids") 
track_indexes = np.load(f"{data_directory}/{ids_file_name}")

print("Remove spectra with many indefinite values", end="\n")

number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
print(f"Indefinite fluxes before drop: {number_indefinite_values}", end="\n")

drop_fraction_spectra = parser.getfloat("processing", "drop_spectra")

keep_spectra_mask = inputting.drop_spectra(
    spectra=spectra, drop_fraction=drop_fraction_spectra
)

spectra = spectra[keep_spectra_mask, :]
specobjids = track_indexes[keep_spectra_mask, 1].reshape(-1, 1)
indexes = np.arange(0, specobjids.size, 1).reshape(-1, 1)

track_indexes = np.hstack((indexes, specobjids))

# update meta data with remaining galaxies
spectra_df = spectra_df.loc[specobjids[:, 0]]
spectra_df.to_csv(f"{data_directory}/drop_{spectra_df_name}")

number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
print(f"Indefinite fluxes after drop: {number_indefinite_values}")

print("Remove wavelegths with many indefinite values", end="\n")

number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
print(f"Indefinite values before drop: {number_indefinite_values}", end="\n")

drop_fraction_waves = parser.getfloat("processing", "drop_waves")

keep_waves_mask = inputting.drop_waves(
    spectra=spectra, drop_fraction=drop_fraction_waves
)

spectra = spectra[:, keep_waves_mask]

print(f"Set new wavelength grid", end="\n")

interpolation_config_file = parser.get("files", "interpolation_config")
interpolation_parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read(interpolation_config_file)

config = ConfigurationFile()

grid_parametes = config.section_to_dictionary(
    interpolation_parser.items("grid")
)

wave = np.linspace(
    grid_parametes["lower"],
    grid_parametes["upper"],
    grid_parametes["number_waves"]
)

wave = wave[keep_waves_mask]

np.save(f"{data_directory}/wave.npy")

print("Inputting indefinite values by the median", end="\n")
nan_median = np.nanmedian(spectra, axis=1)
indefinite_values_mask = ~np.isfinite(spectra)
spectra[indefinite_values_mask] = nan_median
print("Normalize by the median", end="\n")
spectra *= 1 / nan_median.reshape(-1, 1)

np.save(f"{data_directory}/spectra.npy", spectra.astype(np.float32))

print("Save configuration file", end="\n")

with open(f"{data_directory}/{name_config_file}", "w") as configfile:
    parser.write(configfile)

finish_time = time.time()
print(f"Running time: {finish_time - start_time:.2f} [s]")
