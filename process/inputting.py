"""Imputting of missing values in interpolated spectra"""
from configparser import ConfigParser, ExtendedInterpolation
import time

import numpy as np
import pandas as pd

from sdss.process import inputting
from sdss.utils.configfile import ConfigurationFile


start_time = time.time()

parser = ConfigParser(interpolation=ExtendedInterpolation())
name_config_file = "inputting.ini"
parser.read(f"{name_config_file}")

# A load data frame with meta data
data_directory = parser.get("directory", "data")

spectra_df_name = parser.get("files", "spectra_df")
spectra_df = pd.read_csv(
    f"{data_directory}/{spectra_df_name}", index_col="specobjid"
)

# Load interpolated spectra
spectra_file_name = parser.get("files", "spectra")
print(f"{data_directory}/{spectra_file_name}")
spectra = np.load(f"{data_directory}/{spectra_file_name}")
# Load indexes and specobjid of interpolated spectra
ids_file_name = parser.get("files", "ids")
track_indexes = np.load(f"{data_directory}/{ids_file_name}")

variance_file_name = parser.get("files", "variance")
variance_of_spectra = np.load(f"{data_directory}/{variance_file_name}")
#########################################################################
print("Remove spectra with many indefinite values")

drop_fraction_spectra = parser.getfloat("processing", "drop_spectra")

keep_spectra_mask = inputting.drop_spectra(
    spectra=spectra, drop_fraction=drop_fraction_spectra
)

print("Spectra shape", spectra.shape)
spectra = spectra[keep_spectra_mask]
print("Spectra shape", spectra.shape)

specobjids = track_indexes[keep_spectra_mask, 1].reshape(-1, 1)
indexes = np.arange(0, specobjids.size, 1).reshape(-1, 1)
track_indexes = np.hstack((indexes, specobjids))
np.save(f"{data_directory}/ids_inputting.npy", track_indexes)

variance_of_spectra = variance_of_spectra[keep_spectra_mask, :]

# update meta data with remaining galaxies
spectra_df = spectra_df.loc[specobjids[:, 0]]
spectra_df.to_csv(f"{data_directory}/drop_{spectra_df_name}")
#########################################################################
print("Remove wavelegths with many indefinite values")

drop_fraction_waves = parser.getfloat("processing", "drop_waves")

keep_waves_mask = inputting.drop_waves(
    spectra=spectra, drop_fraction=drop_fraction_waves
)

print("Spectra shape", spectra.shape)
spectra = spectra.T
spectra = spectra[keep_waves_mask]
spectra = spectra.T
print("Spectra shape", spectra.shape)
# Save variance of spectra after indefinite values removal
variance_of_spectra = variance_of_spectra[:, keep_waves_mask]
np.save(
    f"{data_directory}/inputting_variance_spectra.npy", variance_of_spectra
)

print("Set new wavelength grid")

interpolation_config_file = parser.get("files", "interpolation_config")
interpolation_parser = ConfigParser(
    interpolation=ExtendedInterpolation()
)
interpolation_parser.read(interpolation_config_file)

config = ConfigurationFile()

grid_parametes = config.section_to_dictionary(
    interpolation_parser.items("grid"), value_separators=[" "]
)

wave = np.linspace(
    grid_parametes["lower"],
    grid_parametes["upper"],
    grid_parametes["number_waves"],
)
print(wave.shape, keep_waves_mask.shape)
wave = wave[keep_waves_mask]

np.save(f"{data_directory}/wave.npy", wave)
#########################################################################
print("Normalize by the median")

median_flux = np.nanmedian(spectra, axis=1)
spectra *= 1 / median_flux.reshape(-1, 1)
#########################################################################
print("Inputting indefinite values by the median")

spectra = inputting.missing_wave_to_median(spectra)
np.save(f"{data_directory}/spectra.npy", spectra.astype(np.float32))
#########################################################################
print("Save configuration file")

with open(
    f"{data_directory}/{name_config_file}", "w", encoding="utf-8"
) as configfile:
    parser.write(configfile)

finish_time = time.time()
print(f"Running time: {finish_time - start_time:.2f} [s]")
