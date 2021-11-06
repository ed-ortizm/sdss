
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
    parser.read("describe.ini")
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
    data_description = spectra_df[["z", "snMedian"]].describe()
    ###########################################################################
    finish_time = time.time()
    print(f"Run time: {finish_time - start_time:.2f}")
