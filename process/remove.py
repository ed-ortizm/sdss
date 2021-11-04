#! /usr/bin/env python3
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

from sdss.raw.data import RawData

###############################################################################
if __name__ == "__main__":
    mp.set_start_method("spawn")

    ti = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("remove.ini")
    meta_data_directory = parser.get("directories", "meta_data")
    ###########################################################################
    keep_df_name = parser.get("files", "keep_df")

    keep_df = pd.read_csv(
        f"{meta_data_directory}/{keep_df_name}", index_col="specobjid"
    )

    ##############################################################
    remove_df_name = parser.get("files", "remove_df")
    remove_df = pd.read_csv(
        f"{meta_data_directory}/{remove_df_name}", index_col="specobjid"
    )
    ##############################################################
    is_in_keep_mask = remove_df.index.isin(keep_df.index)
    remove_mask = ~is_in_keep_mask

    files_remove_df = remove_df.loc[remove_mask]
    ##############################################################
    data_directory = parser.get("directories", "data")
    number_processes = parser.getint("parameters", "number_processes")

    raw = RawData(
        data_directory=data_directory,
        output_directory=data_directory,
        number_processes=number_processes,
    )

    raw.remove_fits_files(files_remove_df)
    ###########################################################################
    tf = time.time()

    print(f"Running time: {tf-ti:.2f} [seg]")
