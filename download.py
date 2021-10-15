#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import time

import numpy as np
import pandas as pd

from src import download

###############################################################################
ti = time.time()
###############################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read("download.ini")
###############################################################################
# SNR sorted data
output_directory = parser.get("directories", "output")
spectra_df_name = parser.get("files", "spectra_df")

spectra_df = pd.read_csv(f"{output_directory}/{spectra_df_name}")

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
################################################################################
tf = time.time()

print(f"Running time: {tf-ti:.2f} [seg]")
