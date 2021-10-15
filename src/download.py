import multiprocessing as mp
import os
import time
import urllib.request

####################################################################
import numpy as np
import pandas as pd

###############################################################################
def init_download_worker(input_counter):
    """
    Initialize worker for download
    PARAMETERS
        counter:
    """
    global counter

    counter = input_counter

    # Update counter to show advance
    # with counter.get_lock():
    #     counter.value += 1
    #     print(f"Download file N: {counter.value}", end="\r")
###############################################################################
class DownloadData:
    def __init__(self,
        spectra_df: "pandas data frame",
        output_directory: "str",
        n_processes: "int"
        ):
        """
        PARAMETERS
            spectra_df: Pandas DataFrame with all the information of
                the sdss galaxies
            output_directory: Path where the data will be downloaded
            n_processes: Number of processes for parallel execution
        """
        self.data_frame = spectra_df

        self._check_directory(output_directory)
        self.output_directory = output_directory

        self.n_processes = n_processes

    ###########################################################################
    def download_files(self):

        print(f"*** Getting {self.data_frame.shape[0]} fits files ****")
        start_time_download = time.time()

        spectra_indexes = self.data_frame.index.values

        counter = mp.Value("i", 0)

        with mp.Pool(
            processes=self.n_processes,
            initializer=init_download_worker,
            initargs=(counter,),
        ) as pool:
            results = pool.map(self._get_file, spectra_indexes)
            number_fail= sum(results)

        finish_time_download = time.time()

        print(f"Done! Finished downloading .fits files...")
        print(f"Failed to download {number_fail} files")
        print(
            f"Download took {finish_time_download-start_time_download:.2f}[s]"
        )

    def _get_file(self, idx_data_frame):

        object = self.data_frame.iloc[idx_data_frame]
        plate, mjd, fiberid, run2d = self._file_identifier(object)

        fname = f"spec-{plate}-{mjd}-{fiberid}.fits"

        SDSSpath = f"sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}"
        folder_path = f"{self.output_directory}/{SDSSpath}"

        url = f"https://data.sdss.org/{SDSSpath}/{fname}"

        if not os.path.exists(folder_path):
            os.makedirs(folder_path, exist_ok=True)

        # Try & Except a failed Download

        try:
            self._retrieve_url(url, folder_path, fname)
            return 0

        except Exception as e:

            print(f"Failed : {url}")

            print(f"{e}")
            return 1

    def _retrieve_url(self, url, folder_path, fname):

        if not (os.path.isfile(f"{folder_path}/{fname}")):

            with counter.get_lock():
                counter.value += 1
                print(f"[{counter.value}] download {fname}", end="\r")
            # print(f"Downloading ", end="\r")

            urllib.request.urlretrieve(url, f"{folder_path}/{fname}")

            file_size = os.path.getsize(f"{folder_path}/{fname}")

            j = 0

            while j < 10 and (file_size < 60000):

                os.remove(f"{folder_path}/{fname}")
                urllib.request.urlretrieve(url, f"{folder_path}/{fname}")
                j += 1
                time.sleep(1)

            file_size = os.path.getsize(f"{folder_path}/{fname}")

            if file_size < 60000:
                print(f"Size of {fname}: {file_size}... Removing file!!")
                os.remove(f"{folder_path}/{fname}")
                raise Exception("Spectra wasn't found")

        else:
            print(f"{fname} already downloaded!!", end="\r")

    def _file_identifier(self, object):

        plate = f"{object['plate']:04}"
        mjd = f"{object['mjd']}"
        fiberid = f"{object['fiberid']:04}"
        run2d = f"{object['run2d']}"

        return plate, mjd, fiberid, run2d

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
