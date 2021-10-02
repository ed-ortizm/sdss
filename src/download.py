import multiprocessing as mp
import os
import time
import urllib.request

####################################################################
import numpy as np
import pandas as pd

################################################################################
class DownloadData:
    def __init__(self, files_data_frame, download_path, n_processes):
        """
        files_data_frame: Pandas DataFrame with all the imformation of the sdss
        galaxies

        download_path: (string) Path where the data will be downloaded
        """
        self.data_frame = files_data_frame
        self.download_path = download_path
        self.n_processes = n_processes

    def get_files(self):

        print(f"*** Getting {len(self.data_frame)} fits files ****")
        start_time_download = time.time()

        if not os.path.exists(self.download_path):
            os.makedirs(self.download_path)

        params = range(len(self.data_frame))

        with mp.Pool(processes=self.n_processes) as pool:
            res = pool.map(self._get_file, params)
            n_failed = sum(res)

        finish_time_download = time.time()

        print(f"Done! Finished downloading .fits files...")
        print(f"Failed to download {n_failed} files")
        print(
            f"Download took {finish_time_download-start_time_download:.2f}[s]"
        )

    def _get_file(self, idx_data_frame):

        object = self.data_frame.iloc[idx_data_frame]
        plate, mjd, fiberid, run2d = self._file_identifier(object)

        fname = f"spec-{plate}-{mjd}-{fiberid}.fits"

        SDSSpath = f"sas/dr16/sdss/spectro/redux/{run2d}/spectra/lite/{plate}"
        folder_path = f"{self.download_path}/{SDSSpath}"

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

            print(f"Downloading {fname}", end="\r")

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
