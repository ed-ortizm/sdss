"""
Take all raw spectra and count the number of fluxes
that are masked because the std of it is larger than
the measured flux
"""
import argparse
from configparser import ConfigParser, ExtendedInterpolation
import time

import numpy as np
import pandas as pd

def count_std_larger_flux(specobjid, data_dir):

    """
    Count number of masked fluxes based on varaince
    """

    wave_flux_ivar = np.load(
        f'{data_dir}/{specobjid}.npy'
    )

    flux = wave_flux_ivar[1, :]
    ivar = wave_flux_ivar[2, :]
    # Get variance of each flux
    ivar[ivar == 0] = np.nan
    variance = 1 / ivar
    variance[np.isnan(ivar)] = np.inf
    variance = np.nan_to_num(variance)

    # mask of: variance > flux [higly uncertain values]
    flux_no_nans = np.nan_to_num(flux)
    large_variance_mask = np.sqrt(variance) > flux_no_nans

    count = large_variance_mask.sum()

    return count

def main():
    """Main"""

    parser = argparse.ArgumentParser(
        description="Find flux/std < 1"
        )

    parser.add_argument(
        "--config",
        type=str,
        default="variance_masking.ini",
        help="Path to config file"
    )

    args = parser.parse_args()

    config_path = args.config

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read(config_path)

    start_time = time.perf_counter()

    # config_handler = ConfigurationFile()

    print("Load data")
    data_dir = parser.get("directories", "data_dir")
    spec_dir = parser.get("directories", "spec_dir")

    spectra_df_name = parser.get("files", "spectra_df")
    spectra_df = pd.read_csv(
        f"{spec_dir}/{spectra_df_name}",
        index_col="specobjid"
    )

    number_spectra = parser.getint(
        "parameters", "number_spectra"
    )

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]

    count_arr = np.empty((spectra_df.shape[0], 2))

    for idx, specobjid in enumerate(spectra_df.index.values):

        count = count_std_larger_flux(specobjid, data_dir)

        count_arr[idx, 0] = specobjid
        count_arr[idx, 1] = count

    np.save(
        f"{spec_dir}/count_std_larger_than_flux.npy",
        count_arr
    )

    finish_time = time.perf_counter()
    print(f"Running time: {finish_time-start_time:.2f}")

if __name__ == "__main__":

    main()
