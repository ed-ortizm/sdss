import multiprocessing as mp
import os

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from sdss.utils.managefiles import FileDirectory

#########################################################################
def init_sample_data_worker(
    input_counter: mp.Value, input_df: pd.DataFrame
) -> None:
    """
    Initialize worker to get sample relevant for the science
    PARAMETERS
        counter: counts the number of the child process
        input_df: pandas dataframe to be accessible to each child
    """
    global counter
    global galaxies_df

    counter = input_counter
    galaxies_df = input_df


#########################################################################
class SampleData(FileDirectory):
    """Use meta data to define a sample for science"""

    def __init__(
        self,
        # data_directory: "str",
        # output_directory: "str",
        # number_processes: "int",
    ):
        """
        PARAMETERS

            # data_directory : sdss raw data's directory
            output_directory : save here .npy files
            number_processes : number of processes to use with mp.Pool

        OUTPUT
            SampleData object
        """

        # super().check_directory(data_directory, exit=True)
        # self.data_directory = data_directory

        # super().check_directory(output_directory, exit=False)
        # self.data_output_directory = output_directory

        # self.number_processes = number_processes

    #####################################################################
    def red_shift(
        self,
        spectra_df: pd.DataFrame,
        lower_bound: float,
        upper_bound: float,
    ) -> pd.Series:

        """
        PARAMETERS

        OUTPUTS
            sample_selection: contains the specobjid to sample from
                original data frame

        """

        z_no_qso_mask = spectra_df["z_noqso"].values != 0.0

        spectra_df.loc[z_no_qso_mask, "z"] = spectra_df.loc[
            z_no_qso_mask, "z_noqso"
        ]

        z = spectra_df["z"]

        sample_selection_mask = (lower_bound < z) * (z < upper_bound)

        return sample_selection_mask

    #####################################################################
    def signal_to_noise(
        self,
        spectra_df: pd.DataFrame,
        lower_bound: float,
        upper_bound: float,
    ) -> np.array:
        """
        PARAMETERS

        OUTPUTS
            sample_selection_mask: contain bools to pick sample from
                original data frame

        """
        signal_to_noise = spectra_df["snMedian"]

        sample_selection_mask = (lower_bound < signal_to_noise) * (
            signal_to_noise < upper_bound
        )

        return sample_selection_mask
