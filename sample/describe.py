#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from sdss.describe import DataDescription
from sdss.process.sample import FileDirectory
from sdss.process.sample import SampleData

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    start_time = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("describe.ini")
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
    str_to_list = lambda x: x.replace(" ", "").split(",")
    ###########################################################################
    describe_parameters = parser.get("parameters", "describe")
    describe_parameters = str_to_list(describe_parameters)

    data_description = spectra_df[describe_parameters].describe()
    ###########################################################################
    # Set parameters for memoir table
    formatter = lambda x: f"{x:.4f}" if x%1 != 0 else f"{x:.0f}"

    save_to = parser.get("latex", "save_to")

    header = parser.get("latex", "header")
    header = str_to_list(header)

    formatters = [formatter, formatter]
    caption = parser.get("latex", "caption")
    label = parser.get("latex", "label")
    bold_rows = parser.getboolean("latex", "bold_rows")
    position = parser.get("latex", "position")
    index = parser.getboolean("latex", "index")

    DataDescription().description_to_latex(
        data_description=data_description,
        save_to=save_to,
        header=header,
        formatters=formatters,
        caption=caption,
        label=label,
        bold_rows=bold_rows,
        position=position,
        index=index,
    )
    print(caption)
    ###########################################################################
    finish_time = time.time()
    print(f"Run time: {finish_time - start_time:.2f}")
