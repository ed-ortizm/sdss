#! /usr/bin/env python3
####################################################################
# Idea is to handle everything from the data frame
####################################################################
from configparser import ConfigParser, ExtendedInterpolation
import time

import numpy as np
import pandas as pd

###############################################################################
start_time = time.time()

parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read("train_sets.ini")
###############################################################################
# A load data frame with meta data
print(f"Load data frame with metadata", end="\n")

meta_data_directory = parser.get("directories", "meta_data")

spectra_df_name = parser.get("files", "spectra_df")
spectra_df = pd.read_csv(
    f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
)

print(f"Load spectra and index array", end="\n")

in_out_directory = parser.get("directories", "output")

data_name = parser.get("files", "spectra")
data = np.load(f"{in_out_directory}/{data_name}")

# this array relates the index of the array with the index in the data frame
index_name = parsed.get("files", "indexes")
index_data = np.load(f"{in_out_directory}/{index_name}")

# Augment data frame with integer position of specobjid in spectra array
spectra_df.loc[index_data[:, 1], "indexArray"] = index_data[:, 0]
###############################################################################
# print(f"Get data frame according to specobjid in spectra array", end="\n")

# specobjid_from_array = index_data[:, 1]
# del index_data

# spectra_df = spectra_df.loc[specobjid_from_array]
###############################################################################
print(f"Save spectra with zWarning", end="\n")

warning_mask = spectra_df["zWarning"] != 0
number_warnings = np.invert(warning_mask).sum()
print(f"Number of warnings: {number_warnings}", end="\n")
index_warning = spectra_df.loc[warning_mask, "indexArray"]
spectra_warnings = np.save("fluxes_with_warnings.npy", data[index_warning])
# get meta data of spectra without flags
spectra_df = spectra_df[np.invert(warning_mask)]
###############################################################################
print(f"Get snMedian bins", end="\n")

spectra_df.sort_values(by=["snMedian"], inplace=True)

number_bins = parser.getint("parameters", "number_bins")

split_number = spectra_df.shape[0] // number_bins
number_remaining_spectra = spectra_df.shape[0] % split_number

data_slices = [split_number * i for i in range(1, number_bins + 1)]

left_slice = 0

for n, right_slice in enumerate(data_slices):

    index_slice = spectra_df["indexArray"].iloc[left_slice:right_slice]

    snr_min = spectra_df.iloc[left_slice]["snMedian"]
    snr_max = spectra_df.iloc[right_slice]["snMedian"]

    array_name = (
        f"bin_{n:02d}_fluxes_"
        f"snr_{snr_min:05.2f}_{snr_max:05.2f}"
    )

    np.save(f"{in_out_directory}/{array_name}.npy", data[index_slice])

    np.save(
        f"{in_out_directory}/{array_name}_shuffle.npy",
        np.random.shuffle(data[index_slice])
    )

    specobjid_slice = spectra_df.iloc[left_slice:right_slice].index
    index_specobjid_slice = np.stack(index_slice, specobjid_slice, axis=1)

    array_name = (
        f"bin_{n:02d}_index_specobjid_"
        f"snr_{snr_min:05.2f}_{snr_max:05.2f}.npy"
    )

    np.save(f"{in_out_directory}/{array_name}", index_specobjid_slice)

    left_slice = right_slice
###############################################################################
# Get remaining slices

if number_remaining_spectra > 1:

    index_slice = spectra_df["indexArray"].iloc[-number_remaining_spectra:]

    snr_min = spectra_df.iloc[-number_remaining_spectra]["snMedian"]
    snr_max = spectra_df.iloc[-1]["snMedian"]

    array_name = (
        f"bin_last_fluxes_"
        f"snr_{snr_min:05.2f}_{snr_max:05.2f}"
    )

    np.save(f"{in_out_directory}/{array_name}.npy", data[index_slice])

    np.save(
        f"{in_out_directory}/{array_name}_shuffle.npy",
        np.random.shuffle(data[index_slice])
    )

    array_name = (
        f"bin_last_index_specobjid_"
        f"snr_{snr_min:05.2f}_{snr_max:05.2f}.npy"
    )

    np.save(f"{in_out_directory}/{array_name}", index_specobjid_slice)
###############################################################################

finish_time = time.time()
print(f"Running time: {finish_time - start_time:.2f} [s]")
