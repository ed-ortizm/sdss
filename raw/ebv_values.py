
"""Spectra processing"""
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
from multiprocessing.sharedctypes import Array
import time

import pandas as pd

from sdss.process import deredspectra

###############################################################################
if __name__ == "__main__":

    mp.set_start_method("spawn")

    start_time = time.time()

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    name_config_file = "ebv_values.ini"
    parser.read(f"{name_config_file}")

    # Load data frame with meta data
    meta_data_directory = parser.get("directories", "meta_data")

    meta_data_name = parser.get("files", "meta_data")
    meta_data = pd.read_csv(
        f"{meta_data_directory}/{meta_data_name}", index_col="specobjid"
    )

    # col 1: specobjid
    # col 2: E(B-V) value
    ebv_values = Array(
        typecode_or_type="d",
        size_or_initializer=2*meta_data.shape[0]
    )

    maps_directory = parser.get("directories", "ebv_maps")
    counter = mp.Value("i", 0)

    number_processes = parser.getint("parameters", "processes")
    specobjids = meta_data.index

    with mp.Pool(
        processes=number_processes,
        initializer=deredspectra.shared_ebv_data,
        initargs=(
            ebv_values,
            meta_data,
            maps_directory,
            counter
        ),
    ) as pool:

        pool.map(deredspectra.ebv_worker, specobjids)

    # save data
    specobjids = ebv_values[:, 0].astype(int)
    meta_data.loc[specobjids, "ebv"] = ebv_values[:, 1]
    meta_data.to_csv(f"{meta_data_directory}/{meta_data_name}")

    finish_time = time.time()
    print(f"Running time: {finish_time - start_time:.2f} [s]")
