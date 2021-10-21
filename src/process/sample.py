import multiprocessing as mp
import os

import astropy.io.fits as pyfits
import numpy as np
import pandas as pd

from src.superclasses import FileDirectory

###############################################################################
def init_sample_data_worker(
    input_counter: "mp.Value", input_df: "pandas dataframe"
) -> "None":
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


###############################################################################
class SampleData(FileDirectory):
    """Use meta data to define a sample for science"""

    def __init__(
        self,
        data_directory: "str",
        output_directory: "str",
        number_processes: "int",
    ):
        """
        PARAMETERS

            data_directory : sdss raw data's directory
            output_directory : save here .npy files
            number_processes : number of processes to use with mp.Pool

        OUTPUT
            GetRawData object
        """

        self.number_processes = number_processes

        super().check_directory(data_directory, exit=True)
        self.data_directory = data_directory

        super().check_directory(output_directory)
        self.data_output_directory = output_directory

    ###########################################################################
    def sample_red_shift(self, spectra_df: "pandas data frame") -> "":
        pass

    ###########################################################################
    def sample_signal_to_noise_ratio(
        self, spectra_df: "pandas data frame"
    ) -> "":

        pass


###############################################################################
