"""
Functionality to handle indefinite values of interpolated spectra
"""

import numpy as np


def missing_wave_to_median(spectra: np.array) -> np.array:
    """
    Replace each missing wavelength with its median value
    throughout the whole array of spectra

    INPUT
    spectra: array of spectra with missing values
        (n_spectra, fluxes)

    OUTPUT
    spectra: array of spectra with missing values replaced
        by its median in each wavelength
    """

    median_flux = np.nanmedian(spectra, axis=0)

    for idx, fluxes in enumerate(spectra.T):

        missing_mask = ~np.isfinite(fluxes)

        fluxes[missing_mask] = median_flux[idx]

    return spectra


def drop_spectra(spectra: np.array, drop_fraction: float) -> np.array:

    """ "
    Drop spectra with a fraction of indefinite values larger
    than the provided 'drop_fraction' argument

    INPUTS
    spectra: array with interplated spectra in a common grid
    drop_fraction: threshold to indicate droping mask

    OUTPUT
    keep_spectra_mask: contains what spectra will be droped to use
        for further processing with specobjid arrays
    """

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra), axis=1)

    number_waves = spectra.shape[1]

    drop_threshold = number_waves * drop_fraction

    keep_spectra_mask = number_indefinite_values < drop_threshold

    return keep_spectra_mask


def drop_waves(spectra: np.array, drop_fraction: float) -> np.array:

    """
    Drop wavelengths in spectra array with a fraction of indefinite
    values larger than the provided 'drop_fraction' argument

    INPUTS
    spectra: array with interplated spectra in a common grid
    drop_fraction: threshold to indicate droping mask

    OUTPUT
    keep_waves_mask: contains what spectra will be droped to use
        for further processing with specobjid arrays
    """

    number_indefinite_values = np.count_nonzero(~np.isfinite(spectra), axis=0)

    number_spectra = spectra.shape[0]

    drop_threshold = number_spectra * drop_fraction

    keep_waves_mask = number_indefinite_values < drop_threshold

    return keep_waves_mask
