#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sdss.describe import DataDescription
from sdss.process.sample import FileDirectory
from sdss.process.sample import SampleData

###############################################################################
start_time = time.time()
###########################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read("describe.ini")
# check files and directories
check = FileDirectory()
###########################################################################
# get relevant data
meta_data_directory = parser.get("directories", "meta_data")

spectra_df_name = parser.get("files", "spectra_df")

spectra_df = pd.read_csv(
    f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
)

number_spectra = parser.getint("parameters", "number_spectra")

if number_spectra != -1:
    spectra_df = spectra_df[:number_spectra]
###############################################################################
str_to_list = lambda x: x.replace(" ", "").split(",")
###############################################################################
generate_stats = parser.getboolean("stats", "generate")

if generate_stats:

    describe_parameters = parser.get("stats", "variables")
    describe_parameters = str_to_list(describe_parameters)

    data_description = spectra_df[describe_parameters].describe()
    ###########################################################################
    # Set parameters for memoir table
    formatter = lambda x: f"{x:.4f}" if x % 1 != 0 else f"{x:.0f}"

    save_to = parser.get("stats", "save_to")

    header = parser.get("stats", "header")
    header = str_to_list(header)

    formatters = [formatter, formatter]
    caption = parser.get("stats", "caption")
    label = parser.get("stats", "label")
    bold_rows = parser.getboolean("stats", "bold_rows")
    position = parser.get("stats", "position")
    index = parser.getboolean("stats", "index")

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
generate_histogram = parser.getboolean("histogram", "generate")

if generate_histogram:

    save_to = parser.get("histogram", "save_to")
    check.check_directory(save_to, exit=False)

    variables = parser.get("histogram", "variables")
    variables = str_to_list(variables)

    number_bins = parser.getint("histogram", "bins")

    for variable in variables:

        spectra_df[variable].hist(bins=number_bins, label=variable)

        location = f"{save_to}/{variable}.pdf"
        plt.savefig(location)
###############################################################################
finish_time = time.time()
print(f"Run time: {finish_time - start_time:.2f}")
