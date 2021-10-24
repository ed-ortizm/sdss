import os
import sys

###############################################################################
class FileDirectory:
    """Handle common operations with files and directories"""

    def __init__(self):
        pass

    ###########################################################################
    def check_directory(
        self, directory: "str", exit: "bool" = False
    ) -> "None":
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
    def file_exists(self, location: "str", exit: "bool" = False) -> "bool":
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
    def remove_file(self, file_location: "str")->"None":

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
    def science_archive_server_url(self):
        pass
    ###########################################################################

    def get_file_location_sas(self, file_row: "pd.row") -> "list":
        """
        PARAMETERS
            file_row:
        OUTPUTS
            return [file_directory, spectrum_name]
                file_directory: location of the spectrum fits file
                spectrum_name: f'spec-{plate}-{mjd}-{fiberid}'
        """


        [plate, mjd, fiberid, run2d] = self._galaxy_identifiers(file_row)

        spectrum_name = f"spec-{plate}-{mjd}-{fiberid}"

        file_directory = (
            f"{self.data_directory}/sas/dr16/sdss/spectro/redux"
            f"/{run2d}/spectra/lite/{plate}"
        )

        return [file_directory, spectrum_name]

    ###########################################################################
    def galaxy_identifiers(self, file_row: "df.row") -> "list":
        """
        PARAMETER
            file_row : pd.row from the object data frame passed to
                the constructor of the class

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
