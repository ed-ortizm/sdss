"""Module to handle common operations with files and directories"""
import os
import sys


class FileDirectory:
    """Handle common operations with files and directories"""

    def __init__(self):
        pass

    @staticmethod
    def check_directory(directory: str, exit_program: bool = False) -> None:
        """
        Check if a directory exists, if not it creates it or
        exits depending on the value of exit
        """

        if os.path.exists(directory) is False:

            if exit_program is True:
                print(f"Directory {directory} NOT FOUND")
                print("Code cannot execute")
                sys.exit()

            os.makedirs(directory)

    @staticmethod
    def file_exists(location: str, exit_program: bool = False) -> bool:
        """
        Check if a location is a file, if not exits depending
        on the value of exit
        """

        file_exists = os.path.isfile(location)

        if file_exists is False:

            file_name = location.split("/")[-1]

            if exit_program is True:
                print(f"File {file_name} NOT FOUND!")
                print("Code cannot execute")
                sys.exit()

            return file_exists

        return file_exists

    @staticmethod
    def remove_file(file_location: str) -> None:

        """
        remove file at file_location
        PARAMETERS
            file_location: e.g. "/home/user/file.text"
        """

        file_name = file_location.split("/")[-1]

        is_file = FileDirectory.file_exists(file_location, exit_program=False)

        if is_file is False:

            print(f"There is no {file_name} at {file_location}")

        else:
            os.remove(file_location)
            print(f"File {file_name} removed!", end="\r")
