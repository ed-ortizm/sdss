"""Module to handle meta data from SDSS spectra"""
import urllib.request

import numpy as np
import pandas as pd


class MetaData:
    """Deal with medata data"""

    def __init__(self):
        pass

    @staticmethod
    def get_z_warning_meaning(warning_flag: int) -> list:

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
    @staticmethod
    def get_sdss_image(
        specobjid: int,
        coordinates: tuple,
        save_to: str,
        image_format: str,
        dimensions: tuple = (0.2, 200, 200),
    ) -> None:

        """
        Download sdss image of the galaxy associated to specobjid

        PARAMETERS
            specobjid:
            RA:
            DEC:
            save_to:
            image_format:
            scale:
            width:
            height:
            image_format:

        """

        sdss_url = (
            f"http://skyserver.sdss.org/dr16/SkyServerWS/ImgCutout/"
            f"getjpeg?TaskName=Skyserver.Explore.Image"
        )

        RA, DEC = coordinates
        coordinates = f"ra={RA}&dec={DEC}"

        scale, width, height = dimensions
        image_dimensions = f"scale={scale}&width={width}&height={height}"

        image_url = f"{sdss_url}&{coordinates}&{image_dimensions}&opt=G"

        urllib.request.urlretrieve(
            image_url, f"{save_to}/{specobjid}.{image_format}"
        )

    ###########################################################################
    @staticmethod
    def download_sdss_spectrum_image(
        specobjid: int, save_to: str, image_format: str
    ) -> None:
        """
        PARAMETERS
            specobjid: unique identification of a galaxy in
                the data frame with the meta data
            save_to: directory location to save the image
            image_format:
        """

        sdss_url = f"http://skyserver.sdss.org/dr16/en/get/SpecById.ashx?id="

        spectrum_url = f"{sdss_url}{specobjid}"

        urllib.request.urlretrieve(
            spectrum_url, f"{save_to}/{specobjid}.{image_format}"
        )

    ###########################################################################
    @staticmethod
    def get_sky_server_url(specobjid: int) -> str:
        """
        PARAMETERS
            specobjid: unique identification of a galaxy in
                the data frame with the meta data

        OUTPUTS
            galaxy_url: url of sdss dr16 in the object explorer
        """

        explorer_url = "http://skyserver.sdss.org/dr16/en/tools/explore"
        galaxy_id = f"summary.aspx?sid={specobjid}&apid="

        galaxy_url = f"{explorer_url}/{galaxy_id}"

        return galaxy_url

    ###########################################################################
    @staticmethod
    def get_file_location_sas(file_row: pd.Series) -> list:
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

        [plate, mjd, fiberid, run2d] = MetaData.galaxy_identifiers(file_row)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"sas/dr16/sdss/spectro/redux" f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectrum_name]

    @staticmethod
    def galaxy_identifiers(file_row: pd.Series) -> list:
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
