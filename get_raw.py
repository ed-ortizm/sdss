from configparser import ConfigParser, ExtendedInterpolation
import numpy as np
import os

import astropy.io.fits as pyfits
import pandas as pd
###############################################################################
def check_file(location: "str", exit: "bool" = False):
    """
    Check if a location exists, if not it creates it or
    exits depending on the value of exit
    """

    file_exists = not os.path.exists(location)

    if not file_exists and exit:

        print(f"File {location} NOT FOUND!")
        print("Code cannot execute")
        sys.exit()

    return file_exists
###############################################################################
parser = ConfigParser(interpolation=ExtendedInterpolation())
parser.read("raw.ini")
# (ReadDataFrame)
raw_frame_location = parser.get("files", "raw_meta_data")
raw_data_frame = pd.read_csv(raw_frame_location)
print(f"Entries in data frame: {raw_data_frame.shape}")
####################################################################
# Use z_noqso when available
raw_data_frame.z = np.where(
                            raw_data_frame.z_noqso.ne(0),
                            raw_data_frame.z_noqso,
                            raw_data_frame.z
                        )

# Remove galaxies with redshift z<=0.01
raw_data_frame = raw_data_frame[raw_data_frame.z > 0.01]
raw_data_frame.index = np.arange(raw_data_frame.shape[0])

print(f"Entries in data frame after z removal: {raw_data_frame.shape}")
####################################################################
number_spectra = parser.getint("parameters", "number_spectra")

if number_spectra != -1:
    raw_data_frame = raw_data_frame[:number_spectra]
####################################################################
data_directory = parser.get("directories", "data")
output_directory = parser.get("directories", "output")
# get raw spectra
classification_data = open(f"{output_directory}/classification.csv", "a")
classification_data.write(f"index,class,subclass\n")

for index in raw_data_frame.index:

    plate = f"{raw_data_frame['plate'].iloc[index]:04}"
    mjd = f"{raw_data_frame['mjd'].iloc[index]}"
    fiberid = f"{raw_data_frame['fiberid'].iloc[index]:04}"
    run2d = f"{raw_data_frame['run2d'].iloc[index]}"

    spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

    print(f"Get data for {spectrum_name} --> {index}", end="\r")

    file_location =(
        f"{data_directory}/sas/dr16/sdss/spectro/redux"
        f"/{run2d}/spectra/lite/{plate}"
        f"/{spectrum_name}.fits"
    )

    if check_file(location=file_location):

        classification_data.write(f"{index},no,no\n")
        continue

    with pyfits.open(file_location) as hdul:

        wave = 10.0 ** (hdul[1].data["loglam"])
        flux = hdul[1].data["flux"]
        classification = hdul[2].data["CLASS"][0]
        sub_class = hdul[2].data["SUBCLASS"][0]

    save_to = f"{output_directory}/rest_frame/{spectrum_name}.npy"
    array_to_save = np.vstack((wave, flux))

    np.save(save_to, array_to_save)

    classification_data.write(f"{index},{classification},{sub_class}\n")

classification_data.close()
