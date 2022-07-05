"""E(B-V) values to correct dust effects on spectra fluxes"""

import multiprocessing as mp
from multiprocessing.sharedctypes import Array

# import numpy as np
import pandas as pd
import sfdmap

from sdss.utils.parallel import to_numpy_array


def get_ebv_value(
    right_ascention: float, declination: float, maps_directory:str
)-> float:
    """
    Compute E(B-V) values from Schlegel, Finkbeiner & Davis (1998)
    dust map FITS files: https://github.com/kbarbary/sfdmap

    INPUTS
    right_ascention: ra of the galaxy
    declination: dec of the galaxy
    maps_directory: location of fits files with E(B-V) maps

    OUTPUT
    ebv_value: value of E(B-V) at (ra, dec) from
        Schlegel, Finkbeiner & Davis (1998)
    """
    ebv_map = sfdmap.SFDMap(maps_directory)
    ebv_value = ebv_map.ebv(right_ascention, declination)

    return ebv_value

def shared_ebv_data(
    ebv_values: Array,
    meta_data: pd.DataFrame,
    maps_directory:str,
    counter: mp.Value
):
    """share data among child processes"""

    global shared_maps_directory
    global ebv_counter
    global shared_ebv_values
    global shared_meta_data

    ebv_counter = counter
    shared_maps_directory = maps_directory

    # first column for specobjid and second for ebv value
    shared_ebv_values = to_numpy_array(
        input_array=ebv_values, array_shape=(meta_data.shape[0], 2)
    )

    shared_meta_data = meta_data

def ebv_worker(specobjid: int) -> None:
    """Obtain E(B-V) values"""

    right_ascention = shared_meta_data.loc[specobjid]["ra"]
    declination = shared_meta_data.loc[specobjid]["dec"]

    with ebv_counter.get_lock():

        ebv_value = get_ebv_value(
            right_ascention, declination, shared_maps_directory
        )

        shared_ebv_values[ebv_counter.value, :] = specobjid, ebv_value

        ebv_counter.value += 1
