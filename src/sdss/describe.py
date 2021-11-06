import pandas as pd
###############################################################################
class DataDescription:
    def __init__(self):
        pass

    def description_to_latex(self,
        data_description: pandas.DataFrame,
        save_to: str,
        header: list,
        formatters: list,
        caption: tuple,
        label: str,
        bold_rows: bool=True,
        position: str="!ht"
        index: bool=True,
    )-> None:
        """
        Converts descriptive statistics of pandas to a latex table
        that uses \usepackage{booktabs} to be properly processed

        PARAMETERS
            data_description:
            save_to:
            header:
            index:
            formatters:
            caption:
            label:
            bold_rows:
            position:

        """

        data_description.to_latex(
            buf=save_to,
            header=header,
            formatters=formatters,
            caption=caption,
            label=label,
            bold_rows=bold_rows,
            position=position,
            index=index
        )
###############################################################################
