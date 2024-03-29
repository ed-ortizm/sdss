"""Extract relevant data from fits files"""
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time

import pandas as pd

from sdss.raw import data

###############################################################################
# spawn creates entirely new processes independent from the parent process
# fork [default] basically just does a minimal cloning, keeping a lot of
# shared elements

# When using spawn you should guard the part that launches
# the job in if __name__ == '__main__':
# set_start_method should also go there

if __name__ == "__main__":
    mp.set_start_method("spawn")

    t0 = time.time()
    ###########################################################################
    parser = ConfigParser(interpolation=ExtendedInterpolation())
    parser.read("raw.ini")
    ###########################################################################
    meta_data_directory = parser.get("directories", "meta_data")

    spectra_df_name = parser.get("files", "spectra_df")
    spectra_df = pd.read_csv(
        f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
    )

    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]
    ###########################################################################
    data_directory = parser.get("directories", "data")
    output_directory = parser.get("directories", "output")
    number_processes = parser.getint("parameters", "number_processes")

    raw_data = data.RawData(
        data_directory=data_directory,
        output_directory=output_directory,
        number_processes=number_processes,
    )
    ###########################################################################
    print("Get raw spectra")
    # raw_data.save_raw_data()
    raw_data.save_raw_data(spectra_df)
    ###########################################################################
    t1 = time.time()
    print(f"Run time: {t1-t0}")
