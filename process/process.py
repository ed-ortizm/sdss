"""Spectra processing"""
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time
from tkinter.tix import INTEGER

import numpy as np
import pandas as pd

from sdss.process import process

if __name__ == "__main__":

    mp.set_start_method("spawn")

    start_time = time.time()

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    name_config_file = "process.ini"
    parser.read(f"{name_config_file}")

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

    # interpolate spectra
    have_to_interpolate = parser.getboolean("parameters", "interpolate")

    if have_to_interpolate is False:

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

    np.save(f"{output_directory}/fluxes.npy", spectra.astype(np.float32))
    np.save(f"{output_directory}/wave.npy", wave)
    ###########################################################################
    # Save configuration file
    with open(f"{output_directory}/{name_config_file}", "w") as configfile:
        parser.write(configfile)
    ###########################################################################
    finish_time = time.time()
    print(f"Running time: {finish_time - start_time:.2f} [s]")

    # Create shared array to store spectra in my desired grid

    # INTERPOLATE
    # (remove means set to nan)
    # (no all raw spectra has the same number of fluxes)
    # remove sky
    # Remove large relative uncertainties (where std>flux)
    # Deredenning spectrum

    # De-redshift spectrum
    # interpolate in common grid
    # filter
    # Removing spectra with many NaNs.

    # Remove fluxes with many nans

    # Normalize spectra by median
    # Impute NaNs by median