"""Interpolate spectra to a common grid in parallel"""
from configparser import ConfigParser, ExtendedInterpolation
from collections import namedtuple
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import time

import numpy as np
import pandas as pd

from sdss.process import interpolate
from sdss.utils.configfile import ConfigurationFile
from sdss.utils.parallel import to_numpy_array


if __name__ == "__main__":

    mp.set_start_method("spawn")

    start_time = time.time()

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    name_config_file = "interpolate.ini"
    parser.read(f"{name_config_file}")

    config_file = ConfigurationFile()

    # A load data frame with meta data
    meta_data_directory = parser.get("directory", "meta_data")

    spectra_df_name = parser.get("files", "spectra_df")
    spectra_df = pd.read_csv(
        f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
    )
    # set number of rows from data frame
    number_spectra = parser.getint("parameters", "number_spectra")

    # for testing purposes
    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]

    number_spectra = spectra_df.shape[0]

    raw_data_directory = parser.get("directory", "raw_spectra")

    # grid paramenters
    grid_parameters = parser.items("grid")
    grid_parameters = config_file.section_to_dictionary(
        grid_parameters, value_separators=[" "]
    )
    # counter to track spectra and link it with specobjid in
    # track_indexes array
    counter = mp.Value("i", 0)
    track_indexes = RawArray("Q", number_spectra * 2)

    # RawArray for spectra
    spectra = RawArray("d", number_spectra * grid_parameters["number_waves"])
    # RawArray for variance_of_spectra
    variance_of_spectra = RawArray(
        "d", number_spectra * grid_parameters["number_waves"]
    )

    shared_arrays_parameters = (
        spectra,
        (number_spectra, grid_parameters["number_waves"]),
        variance_of_spectra,
        (number_spectra, grid_parameters["number_waves"]),
        track_indexes,
        (number_spectra, 2),
    )
    # Set pool of workers
    number_processes = parser.getint("parameters", "processes")

    with mp.Pool(
        processes=number_processes,
        initializer=interpolate.shared_data,
        initargs=(
            counter,
            spectra_df,
            grid_parameters,
            raw_data_directory,
            shared_arrays_parameters,
        ),
    ) as pool:
        # INTERPOLATE
        # (remove means set to nan)
        # (no all raw spectra has the same number of fluxes)
        # remove sky
        # Remove large relative uncertainties (where std>flux)
        # Deredenning spectrum

        # De-redshift spectrum
        # interpolate in common grid

        pool.map(interpolate.worker_interpolation, spectra_df.index)

    output_directory = parser.get("directory", "output")

    spectra = to_numpy_array(spectra, shared_arrays_parameters[1])

    np.save(f"{output_directory}/interpolated_spectra.npy", spectra)

    variance_of_spectra = to_numpy_array(
        variance_of_spectra, shared_arrays_parameters[3]
    )

    np.save(
        f"{output_directory}/interpolated_variance_spectra.npy",
        variance_of_spectra,
    )

    track_indexes = to_numpy_array(track_indexes, shared_arrays_parameters[5])

    np.save(f"{output_directory}/ids_interpolation.npy", track_indexes)

    with open(f"{output_directory}/{name_config_file}", "w") as configfile:
        parser.write(configfile)

    finish_time = time.time()
    print(f"Running time: {finish_time - start_time:.2f} [s]")
