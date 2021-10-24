#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from src import download

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    ti = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("remove.ini")
    ###########################################################################
    data_directory = parser.get("directories", "data")
    galaxy_df_name = parser.get("files", "spectra_df")

    galaxy_df = pd.read_csv(f"{data_directory}/{galaxy_df_name}")

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]
    ##############################################################
    # Data Download
    number_processes = parser.getint("parameters", "number_processes")

    download_spectra = download.DownloadData(
        spectra_df=spectra_df,
        output_directory=output_directory,
        n_processes=number_processes,
    )

    download_spectra.download_files()
    ###########################################################################
    tf = time.time()

    print(f"Running time: {tf-ti:.2f} [seg]")
