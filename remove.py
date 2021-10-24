#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from src.raw.data import GetRawData

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    ti = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("remove.ini")
    ###########################################################################
    data_directory = parser.get("directories", "data")
    files_df_name = parser.get("files", "files_df")

    files_df = pd.read_csv(
                            f"{data_directory}/{files_df_name}",
                            index_col="specobjid"
                        )

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        files_df = files_df[:number_spectra]
    ##############################################################
    # Data Download
    number_processes = parser.getint("parameters", "number_processes")

    raw = GetRawData(
        data_directory=data_directory,
        output_directory=data_directory,
        number_processes=number_processes,
    )

    raw.remove_fits_files(files_df)
    ###########################################################################
    tf = time.time()

    print(f"Running time: {tf-ti:.2f} [seg]")
