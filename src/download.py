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
    def download_files(self)->"None":
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

        download_time = finish_time_download-start_time_download
        print(f"Download took {download_time:.2f}[s]")

    ###########################################################################

    def _get_file(self, index_spectrum: "int")-> "int":
        """
        Retrieve a single file from this science archive server

        PARAMETERS
            index_spectrum: index of the spectrum in spectra data frame

        RETURN
            0 if download is succesful, 1 otherwise
        """

        df_row_spectrum = self.spectra_df.iloc[index_spectrum]
        [plate, mjd, fiberid, run2d] = self._file_identifier(df_row_spectrum)

        file_name = f"spec-{plate}-{mjd}-{fiberid}.fits"

        # sas: science archive server
        sas_location = (
            f"sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}"
        )

        save_to = f"{self.output_directory}/{sas_location}"

        url = f"https://data.sdss.org/{sas_location}/{file_name}"

        sel._check_directory(save_to, exit=False)

        # Try & Except a failed Download

        try:
            self._retrieve_url(url, save_to, file_name)
            return 0

        except Exception as e:

            print(f"Failed : {url}")

            print(f"{e}")
            return 1

    ###########################################################################
    def _retrieve_url(self, url, save_to, file_name):

        if not (os.path.isfile(f"{save_to}/{file_name}")):

            with counter.get_lock():
                counter.value += 1
                print(f"[{counter.value}] download {file_name}", end="\r")
            # print(f"Downloading ", end="\r")

            urllib.request.urlretrieve(url, f"{save_to}/{file_name}")

            file_size = os.path.getsize(f"{save_to}/{file_name}")

            j = 0

            while j < 10 and (file_size < 60000):

                os.remove(f"{save_to}/{file_name}")
                urllib.request.urlretrieve(url, f"{save_to}/{file_name}")
                j += 1
                time.sleep(1)

            file_size = os.path.getsize(f"{save_to}/{file_name}")

            if file_size < 60000:
                print(f"Size of {file_name}: {file_size}... Removing file!!")
                os.remove(f"{save_to}/{file_name}")
                raise Exception("Spectra wasn't found")

        else:
            print(f"{file_name} already downloaded!!", end="\r")

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
