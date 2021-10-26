#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from src.process.sample import SampleData

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    start_time = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("sample.ini")

    ###########################################################################
    # get relevant data
    data_directory = parser.get("directories", "data")

    spectra_df_name = parser.get("files", "spectra_df")

    spectra_df = pd.read_csv(
        f"{data_directory}/{spectra_df_name}",
        usecols=["specobjid", "z", "z_noqso"],
    )

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]
    ###########################################################################
    # Set sample class
    output_directory = parser.get("directories", "output")
    sample = SampleData(data_directory, output_directory)
    ###########################################################################
    # Set sample parameters
    lower_bound = parser.getfloat("redshift", "lower_bound")
    upper_bound = parser.getfloat("redshift", "upper_bound")

    # sample
    z_selection_mask = sample.red_shift_sampling(
        spectra_df=spectra_df, lower_bound=lower_bound, upper_bound=upper_bound
    )
    print(np.count_nonzero(z_selection_mask))
    ###########################################################################
    finish_time = time.time()
    print(f"Run time: {finish_time - start_time:.2f}")
