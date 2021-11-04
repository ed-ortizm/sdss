#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from sdss import download

###############################################################################
# spawn creates entirely new processes independent from the parent process
# fork [default] basically just does a minimal cloning, keeping a lot of
# shared elements

# When using spawn you should guard the part that launches
# the job in if __name__ == '__main__':
# set_start_method should also go there

if __name__ == "__main__":
    mp.set_start_method("spawn")

    ti = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("download.ini")
    ###########################################################################
    meta_data_directory = parser.get("directories", "meta_data")
    spectra_df_name = parser.get("files", "spectra_df")

    spectra_df = pd.read_csv(f"{meta_data_directory}/{spectra_df_name}")

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]
    ##############################################################
    # Data Download
    output_directory = parser.get("directories", "output")
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
