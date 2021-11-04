#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from sdss.process.sample import FileDirectory
from sdss.process.sample import SampleData

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    start_time = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("sample.ini")
    # Check files and directory
    check = FileDirectory()
    ###########################################################################
    # get relevant data
    meta_data_directory = parser.get("directories", "meta_data")

    spectra_df_name = parser.get("files", "spectra_df")

    spectra_df = pd.read_csv(
        f"{meta_data_directory}/{spectra_df_name}",
        index_col="specobjid",
    )

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]
    ###########################################################################
    # Set sample class
    sample = SampleData(
    )
    ###########################################################################
    # Set sample parameters for snr
    z_lower_bound = parser.getfloat("redshift", "lower_bound")
    z_upper_bound = parser.getfloat("redshift", "upper_bound")

    # sample
    z_selection_mask = sample.red_shift(
        spectra_df=spectra_df,
        lower_bound=z_lower_bound,
        upper_bound=z_upper_bound
    )
    ###########################################################################
    # Set sample parameters for snr
    snr_lower_bound = parser.getfloat("signal_to_noise", "lower_bound")
    snr_upper_bound = parser.getfloat("signal_to_noise", "upper_bound")

    # sample
    signal_to_noise_selection_mask = sample.signal_to_noise(
        spectra_df=spectra_df,
        lower_bound=snr_lower_bound,
        upper_bound=snr_upper_bound
    )
    ###########################################################################
    # Save new data frame
    selection_mask = z_selection_mask * signal_to_noise_selection_mask

    z_name = (
        f"{str(z_lower_bound).replace('.', '_')}"
        f"_z_"
        f"{str(z_upper_bound).replace('.', '_')}"
    )

    signal_to_noise_name = (
        f"{str(snr_lower_bound).replace('.', '_')}"
        f"_snr_"
        f"{str(snr_upper_bound).replace('.', '_')}"
    )

    output_directory = parser.get("directories", "output")

    output_directory = f"{output_directory}/{z_name}_{signal_to_noise_name}"
    check.check_directory(output_directory, exit=False)

    spectra_df_name = f"{z_name}_{signal_to_noise_name}.csv.gz"

    spectra_df.loc[selection_mask].to_csv(
        f"{output_directory}/{spectra_df_name}",
        index=True
        )
    ###########################################################################
    finish_time = time.time()
    print(f"Run time: {finish_time - start_time:.2f}")
