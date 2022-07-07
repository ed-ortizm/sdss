"""Spectra processing"""
from configparser import ConfigParser, ExtendedInterpolation
import multiprocessing as mp
import time
from tkinter.tix import INTEGER

import numpy as np
import pandas as pd

from sdss.process import process

if __name__ == "__main__":

    mp.set_start_method("spawn")

    start_time = time.time()

    parser = ConfigParser(interpolation=ExtendedInterpolation())
    name_config_file = "process.ini"
    parser.read(f"{name_config_file}")

    # A load data frame with meta data
    meta_data_directory = parser.get("directories", "meta_data")

    spectra_df_name = parser.get("files", "spectra_df")
    spectra_df = pd.read_csv(
        f"{meta_data_directory}/{spectra_df_name}", index_col="specobjid"
    )
    # set number of rows from data frame
    number_spectra = parser.getint("parameters", "number_spectra")

    if number_spectra != -1:
        spectra_df = spectra_df[:number_spectra]

    raw_data_directory = parser.get("directories", "raw_spectra")

    # Create shared array to store spectra in my desired grid

    #INTERPOLATE
    # (remove means set to nan)
    # (no all raw spectra has the same number of fluxes)
        # remove sky
        # Remove large relative uncertainties (where std>flux)
        # Deredenning spectrum

    # De-redshift spectrum
    # interpolate in common grid
    def interpolate(self, spectra_df: pd.DataFrame) -> np.array:
        """
        Interpolate rest frame spectra from data directory according to
        wave master  and save it to output directory
        OUTPUT
            fluxes: contains interpolate fluxes
        """

        number_spectra = spectra_df.shape[0]
        print(f"Interpolate {number_spectra} spectra")
        number_waves = self.grid.size

        counter = mp.Value("i", 0)

        fluxes = RawArray(ctypes.c_double, number_spectra * number_waves)
        fluxes_shape = (number_spectra, number_waves)

        # array with counter, specobjid and 0 or 1 if processing is
        # successful columns
        track_indexes = RawArray(ctypes.c_uint64, number_spectra * 3)
        track_indexes_shape = (number_spectra, 3)

        spectra_indexes = spectra_df.index.to_numpy()

        with mp.Pool(
            processes=self.number_processes,
            initializer=init_process_worker,
            initargs=(
                counter,
                spectra_df,
                fluxes,
                fluxes_shape,
                track_indexes,
                track_indexes_shape,
            ),
        ) as pool:

            results = pool.map(self._interpolate, spectra_indexes)

        track_indexes = to_numpy_array(track_indexes, track_indexes_shape)
        success_interpolation_mask = ~track_indexes[:, 2].astype(np.bool)
        track_indexes = track_indexes[success_interpolation_mask]
        save_to = f"{self.output_directory}/indexes_interpolate.npy"
        np.save(save_to, track_indexes)

        fluxes = to_numpy_array(fluxes, fluxes_shape)

        fluxes = fluxes[success_interpolation_mask]

        save_to = f"{self.output_directory}/interpolate.npy"
        np.save(save_to, fluxes)

        return fluxes
 
    # filter
    # Removing spectra with many NaNs.

    # Remove fluxes with many nans

    # Normalize spectra by median
    # Impute NaNs by median






    output_directory = parser.get("directories", "output")

    grid_parameters = dict(parser.items("grid"))
    number_processes = parser.getint("parameters", "processes")

    data_process = process.DataProcess(
        raw_data_directory=raw_data_directory,
        output_directory=output_directory,
        grid_parameters=grid_parameters,
        number_processes=number_processes,
    )

    # interpolate spectra
    have_to_interpolate = parser.getboolean("parameters", "interpolate")

    if have_to_interpolate is False:

        spectra_interpolate = parser.get("files", "interpolate")
        spectra = np.load(f"{output_directory}/{spectra_interpolate}")

    else:

        spectra = data_process.interpolate(spectra_df=spectra_df)

    ###########################################################################
    print(f"Handle indefinite values")

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
    print(f"Indefinite fluxes before drop: {number_indefinite_values}")

    drop_fraction = parser.getfloat("parameters", "drop")

    spectra, wave = data_process.drop_indefinite_values(
        spectra=spectra, drop=drop_fraction
    )

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra))
    print(f"Indefinite fluxes after drop: {number_indefinite_values}")

    print(f"Replace missing fluxes and normalize")

    spectra = data_process.replace_missing_fluxes_and_normalize_by_median(
        spectra
    )

    print(f"Save data")

    np.save(f"{output_directory}/fluxes.npy", spectra.astype(np.float32))
    np.save(f"{output_directory}/wave.npy", wave)
    ###########################################################################
    # Save configuration file
    with open(f"{output_directory}/{name_config_file}", "w") as configfile:
        parser.write(configfile)
    ###########################################################################
    finish_time = time.time()
    print(f"Running time: {finish_time - start_time:.2f} [s]")
