import os
import sys
import urllib.request

import numpy as np
import pandas as pd

###############################################################################
class ConfigurationFile:
    """Manage common operation with configuration files.
    For instance:

        [parameters]
        metric = mse, lp, mad

    when reading metric I can inmediately get it as a list if looping
    is needed in addition, I can get a whole section in the configuration
    file as a dictionary and do the necessary transformation to each key,
    value pair

    """

    def __init__(self):
        pass

    ###########################################################################
    def section_to_dictionary(self,
        section_items:tuple,
        value_separators:list,
        # multiple_lines: bool,
        # comma_values:bool,
    ) -> dict:
        """
        Converts a section in the configuration file to a dictionary.
        WARNING: all values are strings.
        If there is a variable with multiple lines, the values of the
        associated key in the dictionary would be a list

        PARAMETERS
            section_items: items in a section of the configuration file

        OUTPUTS
            section_dictionary: items transformed
        """
        section_dictionary = dict(section_items)

        for key, value in section_dictionary.items():

            value = str(value)
            ###################################################################
            if "\n" in value:

                value = value.split("\n")
                section_dictionary[key] = value

            ###################################################################
            for separator in value_separators:

                if separator in value:

                    value = value.split(separator)
                    section_dictionary[key] = value
            ###################################################################

        section_dictionary = self._transform_values_in_dictionary(
                section_dictionary
            )

        return section_dictionary
    ###########################################################################
    def entry_to_list(self, entry: str, entry_type: type, separator: str)-> list:
        """

        PARAMETERS

            entry: a coma separated string
                architecture: 100, 50, 5, 50, 100

            entry_type: either str, float, int or bool
            separator: separator to use when spliting entry string

        OUTPUTS
            entry_list: list of elements in entry with the type
                100 separator 50 separator 5... --> [100, 50, 5, ...]
        """
        entry = entry.strip().split(separator)
        entry = [entry_type(value) for value in entry]

        return entry
    ###########################################################################
    def _transform_values_in_dictionary(self, dictionary: dict):

        for key, value in dictionary.items():

            is_list = type(value) is list
            value = self._transform_values(value, is_list)

            dictionary[key] = value

        return dictionary
    ###########################################################################
    def _transform_values(self, items: str, is_list: bool=False):

        if is_list is True:

            new_items = []

            for string in items:

                value = self._get_value_from_string(string)
                new_items.append(value)

            return new_items

        return self._get_value_from_string(items)

    ###########################################################################
    def _get_value_from_string(self, string: str):

        """
        Get value from string variable, could be: bool, str, int or float
        """
        string = string.strip()
        #######################################################################
        if (string == "True") or (string == "False"):

            return string == "True"

        #######################################################################
        if (string.isalpha() is True) | ("_" in string) is True  :

            return string
        #######################################################################
        numeric_value = float(string)
        if numeric_value.is_integer():

            return int(numeric_value)

        return numeric_value
    ###########################################################################

    ###########################################################################
###############################################################################
class FileDirectory:
    """Handle common operations with files and directories"""

    def __init__(self):
        pass

    ###########################################################################
    def check_directory(self, directory: str, exit: bool = False) -> "None":
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
    def file_exists(self, location: str, exit: bool = False) -> bool:
        """
        Check if a location is a file, if not exits depending
        on the value of exit
        """

        file_exists = os.path.isfile(location)

        if not file_exists:

            file_name = location.split("/")[-1]

            if exit:
                print(f"File {file_name} NOT FOUND!")
                print("Code cannot execute")
                sys.exit()

            return file_exists

        return file_exists

    ###########################################################################
    def remove_file(self, file_location: str) -> "None":

        """
        remove file at file_location
        PARAMETERS
            file_location: e.g. "/home/user/file.text"
        """

        file_name = file_location.split("/")[-1]

        if not self.file_exists(file_location, exit=False):

            print(file_location)
            print(f"There is no {file_name}!", end="\r")

        else:
            os.remove(file_location)
            print(f"File {file_name} removed!", end="\r")


###############################################################################
class MetaData:
    """Deal with medata data """

    def __init__(self):
        pass

    ###########################################################################
    def get_z_warning_meaning(self, warning_flag: int) -> list:

        """
        Takes the interger representing the z warning and converts it to
            the string indicating the meaning of the warning:

        PARAMETERS
            warning_flag: value between 0 and 9
        OUTPUT
            warning_meaning: list of strings  with the meaning of warning
        """

        map_warning = {
            "0": "good",
            "1": "LITTLE_COVERAGE",
            "2": "SMALL_DELTA_CHI2",
            "3": "NEGATIVE_MODEL",
            "4": "MANY_OUTLIERS",
            "5": "Z_FITLIMIT",
            "6": "NEGATIVE_EMISSION",
            "7": "UNPLUGGED",
            "8": "BAD_TARGET",
            "9": "NODATA",
        }

        if warning_flag == 0:
            return []

        warning_flag = list(np.binary_repr(warning_flag, 9))

        warnings = [
            int(val) * (8 - idx) for idx, val in enumerate(warning_flag)
        ]

        warning_meaning = [map_warning[str(i)] for i in warnings if i != 0]

        return warning_meaning

    ###########################################################################
    def get_sdss_image(
        self,
        galaxy_specobjid: int,
        RA: float,
        DEC: float,
        save_to: str,
        format: str,
        scale: float = 0.2,
        width: int = 200,
        height: int = 200,
    ) -> None:

        """
        Download sdss image of the galaxy associated to galaxy_specobjid

        PARAMETERS
            galaxy_specobjid:
            RA:
            DEC:
            save_to:
            format:
            scale:
            width:
            height:
            format:

        """

        sdss_url = (
            f"http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/"
            f"getjpeg?TaskName=Skyserver.Explore.Image"
        )

        coordinates = f"ra={RA}&dec={DEC}"

        image_dimensions = f"scale={scale}&width={width}&height={height}"

        image_url = f"{sdss_url}&{coordinates}&{image_dimensions}&opt=G"

        urllib.request.urlretrieve(
            image_url, f"{save_to}/{galaxy_specobjid}.{format}"
        )

    ###########################################################################
    def download_sdss_spectrum_image(
        self, galaxy_specobjid: int, save_to: str, format: str
    ) -> None:
        """
        PARAMETERS
            galaxy_specobjid: unique identification of a galaxy in
                the data frame with the meta data
            save_to: directory location to save the image
            format:
        """

        sdss_url = f"http://skyserver.sdss.org/dr16/en/get/SpecById.ashx?id="

        spectrum_url = f"{sdss_url}{galaxy_specobjid}"

        urllib.request.urlretrieve(
            spectrum_url, f"{save_to}/{galaxy_specobjid}.{format}"
        )

    ###########################################################################
    def get_sky_server_url(self, galaxy_specobjid: int) -> str:
        """
        PARAMETERS
            galaxy_specobjid: unique identification of a galaxy in
                the data frame with the meta data

        OUTPUTS
            galaxy_url: url of sdss dr16 in the object explorer
        """

        explorer_url = "http://skyserver.sdss.org/dr16/en/tools/explore"
        galaxy_id = f"summary.aspx?sid={galaxy_specobjid}&apid="

        galaxy_url = f"{explorer_url}/{galaxy_id}"

        return galaxy_url

    ###########################################################################
    def get_file_location_sas(self, file_row: pd.Series) -> list:
        """
        PARAMETERS
            file_row: contains at least the columns
                plate, mjd, fiberid, run2d
        OUTPUTS
            return [file_directory, spectrum_name]
                file_directory: location of the spectrum fits file
                in sas directory
                spectrum_name: f'spec-{plate}-{mjd}-{fiberid}'
        """

        [plate, mjd, fiberid, run2d] = self.galaxy_identifiers(file_row)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"sas/dr16/sdss/spectro/redux" f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectrum_name]

    ###########################################################################
    def get_spectrum_name(self, file_row: pd.Series) -> str:

        [plate, mjd, fiberid, run2d] = self.galaxy_identifiers(file_row)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        pass

    ###########################################################################
    def galaxy_identifiers(self, file_row: "df.row") -> list:
        """
        PARAMETER
            file_row: contains at least the columns
                plate, mjd, fiberid, run2d
        OUTPUT
            return [plate, mjd, fiberid, run2d]
                plate: self explanatory
                mjd: date
                fiberid: self explanatory
                run2d: PENDING

        """

        plate = f"{file_row['plate']:04}"
        mjd = f"{file_row['mjd']}"
        fiberid = f"{file_row['fiberid']:04}"
        run2d = f"{file_row['run2d']}"

        return [plate, mjd, fiberid, run2d]


###############################################################################
