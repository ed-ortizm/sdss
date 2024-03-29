import multiprocessing as mp
import os
import time
import urllib.request

####################################################################
import numpy as np
import pandas as pd

###############################################################################
def init_download_worker(input_counter: "mp.Value") -> "None":
    """
    Initialize worker for download
    PARAMETERS
        counter: counts the number of the child process
    """
    global counter

    counter = input_counter


###############################################################################
class DownloadData:
    def __init__(
        self,
        spectra_df: "pandas data frame",
        output_directory: "str",
        n_processes: "int",
    ) -> "None":
        """
        PARAMETERS
            spectra_df: Pandas DataFrame with all the information of
                the sdss galaxies. Columns of the data frame are:
                [
                    'specobjid', 'mjd', 'plate', 'fiberid', 'run2d',
                    'ra', 'dec', 'z', 'zErr', 'zWarning', 'class',
                    'subClass', 'z_noqso', 'zErr_noqso', 'zWarning_noqso',
                    'snMedian'
                ]
            output_directory: location where the data will be downloaded
            n_processes: number of processes for parallel execution
        """
        self.spectra_df = spectra_df

        self._check_directory(output_directory)
        self.output_directory = output_directory

        self.n_processes = n_processes

    ###########################################################################
    def download_files(self) -> "None":

        """Download spectra from the science archive server"""

        print(f"Start download of {self.spectra_df.shape[0]} fits files...")

        start_time_download = time.time()

        spectra_indexes = self.spectra_df.index.values

        counter = mp.Value("i", 0)

        with mp.Pool(
            processes=self.n_processes,
            initializer=init_download_worker,
            initargs=(counter,),
        ) as pool:

            results = pool.map(self._get_file, spectra_indexes)
            number_fail = sum(results)

        finish_time_download = time.time()

        print(f"Finish download...")
        print(f"Fail to download {number_fail} files")

        download_time = finish_time_download - start_time_download
        print(f"Download took {download_time:.2f}[s]")
        #######################################################################

    ###########################################################################

    def _get_file(self, index_spectrum: "int") -> "int":
        """
        Retrieve a single file from the science archive server

        PARAMETERS
            index_spectrum: index of the spectrum in spectra data frame

        RETURN
            0 if download is succesful, 1 otherwise
        """

        df_row_spectrum = self.spectra_df.iloc[index_spectrum]

        [plate, mjd, fiberid, run2d] = self._file_identifier(df_row_spectrum)

        file_name = f"spec-{plate}-{mjd}-{fiberid}.fits"

        # Try & Except a failed Download

        try:
            self._query_file(file_name, run2d, plate)

            return 0

        except Exception as e:

            print(f"Failed : {file_name}. run2d:{run2d}")
            print(f"{e}")

            return 1

    ###########################################################################
    def _query_file(
        self, file_name: "str", run2d: "str", plate: "str"
    ) -> "None":

        # sas: science archive server
        sas_location = (
            f"sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}"
        )

        save_to = f"{self.output_directory}/{sas_location}"

        file_url = f"https://data.sdss.org/{sas_location}/{file_name}"

        self._check_directory(save_to, exit=False)

        file_exists = self._file_exits(
            file_location=f"{save_to}/{file_name}", exit=False
        )

        if not file_exists:

            with counter.get_lock():
                counter.value += 1
                print(f"[{counter.value}] Download {file_name}", end="\r")

            urllib.request.urlretrieve(file_url, f"{save_to}/{file_name}")

            file_size = os.path.getsize(f"{save_to}/{file_name}")

            self._retry_download_if_small_size(
                file_size, file_name, save_to, file_url
            )

        else:
            print(f"{file_name} already downloaded!!")

    ###########################################################################
    def _retry_download_if_small_size(
        self,
        file_size: "float",
        file_name: "str",
        save_to: "str",
        file_url: "str",
    ) -> "None || exception":

        """
        Check the size of downloaded file. If the size is smaller than
        60 Kbs, it trys to download it again (at least 10 times), otherwise
        it will raise an exception
        """

        j = 0

        while j < 10 and (file_size < 60000):

            os.remove(f"{save_to}/{file_name}")
            urllib.request.urlretrieve(file_url, f"{save_to}/{file_name}")
            j += 1
            time.sleep(1)

        file_size = os.path.getsize(f"{save_to}/{file_name}")

        if file_size < 60000:
            print(f"Size of {file_name}: {file_size}... Removing file!!")
            os.remove(f"{save_to}/{file_name}")
            raise Exception("Spectra wasn't found")

    ###########################################################################
    def _file_identifier(self, df_row_spectrum):

        plate = f"{df_row_spectrum['plate']:04}"
        mjd = f"{df_row_spectrum['mjd']}"
        fiberid = f"{df_row_spectrum['fiberid']:04}"
        run2d = f"{df_row_spectrum['run2d']}"

        return [plate, mjd, fiberid, run2d]

    ###########################################################################
    def _check_directory(self, directory: "str", exit: "bool" = False):
        """
        Check if a directory exists, if not it creates it or
        exits depending on the value of exit
        """

        if not os.path.exists(directory):

            if exit:
                print(f"Directory {diretory} NOT FOUND")
                print("Code cannot execute")
                sys.exit()

            os.makedirs(directory)

    ###########################################################################
    def _file_exits(self, file_location: "str", exit: "bool") -> "bool":
        """
        Check if a file exists and depending on exit parameter, it exits
        the program because the file is necessary for computations.

        PARAMETERS

            file_location: file location with extension
                example: /home/john/data/data.txt

            exit: if True, it exits the program
        """

        if not os.path.isfile(file_location):

            if exit:
                print(f"NOT FOUND: {file_location}")
                print(f"Program cannot execute with out this file")
                sys.exit()

            return False

        return True
