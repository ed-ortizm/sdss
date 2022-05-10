"""Module to handle common operations with configuration files"""


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
    def section_to_dictionary(
        self, section_items: tuple, value_separators: list
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

    @staticmethod
    def entry_to_list(entry: str, entry_type: type, separator: str) -> list:
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

            is_list = isinstance(value, list)
            value = self._transform_values(value, is_list)

            dictionary[key] = value

        return dictionary

    ###########################################################################
    def _transform_values(self, items: str, is_list: bool = False):

        if is_list is True:

            new_items = []

            for string in items:

                value = self._get_value_from_string(string)
                new_items.append(value)

            return new_items

        return self._get_value_from_string(items)

    ###########################################################################
    @staticmethod
    def _get_value_from_string(string: str):

        """
        Get value from string variable, could be: bool, str, int or float
        """
        string = string.strip()
        #######################################################################
        if string in ("True", "False"):

            return string == "True"

        #######################################################################
        if (string.isalpha() is True) | ("_" in string) is True:

            return string
        #######################################################################
        numeric_value = float(string)
        if numeric_value.is_integer():

            return int(numeric_value)

        return numeric_value
