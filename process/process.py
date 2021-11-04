#! /usr/bin/env python3
####################################################################
# Idea is to handle everything from the data frame
####################################################################
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import os
import time

import numpy as np
import pandas as pd

from sdss.process import process

###############################################################################
if __name__ == "__main__":

    mp.set_start_method("spawn")

    start_time = time.time()

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("process.ini")
    ###########################################################################
    raw_data_directory = parser.get("directories", "raw_spectra")
    output_directory = parser.get("directories", "output")

    grid_parameters = dict(parser.items("grid"))
    number_processes = parser.getint("parameters", "processes")

    data_process = process.DataProcess(
        raw_data_directory=raw_data_directory,
        output_directory=output_directory,
        grid_parameters=grid_parameters,
        number_processes=number_processes,
    )
    ###########################################################################
    # A load data frame with meta data
    meta_data_directory = parser.get("directories", "meta_data")

    spectra_df_name = parser.get("files", "spectra_df")
    spectra_df = pd.read_csv(
        f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
    )

    # set number of rows from data frame
    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]

    ###########################################################################
    # interpolate spectra
    have_to_interpolate = parser.getboolean("parameters", "interpolate")

    if not have_to_interpolate:

        spectra_interpolate = parser.get("files", "interpolate")
        spectra = np.load(f"{output_directory}/{spectra_interpolate}")

    else:

        spectra = data_process.interpolate(spectra_df=spectra_df)

    ###########################################################################
    print(f"Handle indefinite values")

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
    print(f"Indefinite fluxes before drop: {number_indefinite_values}")

    drop_fraction = parser.getfloat("parameters", "drop")

    spectra, wave = data_process.drop_indefinite_values(
        spectra=spectra, drop=drop_fraction
    )

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
    print(f"Indefinite fluxes after drop: {number_indefinite_values}")

    print(f"Replace missing fluxes and normalize")

    spectra = data_process.replace_missing_fluxes_and_normalize_by_median(
        spectra
    )

    print(f"Save data")

    np.save(f"{output_directory}/fluxes.npy", spectra)
    np.save(f"{output_directory}/wave.npy", wave)
    ###########################################################################
    finish_time = time.time()
    print(f"Running time: {finish_time - start_time:.2f} [s]")
