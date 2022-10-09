"""Handle indefinite values in spectra"""
import numpy as np

#########################################################################
def replace_missing_fluxes_and_normalize_by_median(
    self, spectra: np.array
) -> np.array:
    """ """

    missing_values_mask = ~np.isfinite(spectra)

    nan_median = np.nanmedian(spectra, axis=1)
    null_median_mask = nan_median == 0.0
    nan_median[null_median_mask] = 1.0

    spectra *= 1 / nan_median.reshape((-1, 1))

    spectra[missing_values_mask] = 1.0

    # zero for spectra where nanmedian is null
    spectra[null_median_mask, :] = 0

    return spectra


#########################################################################
def drop_indefinite_values(self, spectra: np.array, drop: float = 0.1) -> list:
    """
    Drops indefinite values per wavelength according to
    the drop fraction specified by the drop parameter
    PARAMETERS
        spectra: contains interpolated fluxes in a common grid
        drop: fraction of indefinite values to discard at a given
            wavelength, e.g. drop = 0.1 will discard a wavelength
            where more than 10% of fluxes are indefinite values
    OUTPUTS
        [spectra, wave]:
            spectra: fluxes without indefinite values
            wave: common grid to all spectra update according
                to wavelengths dropped
    """
    print(f"Shape before discard indefinite values: {spectra.shape}")

    number_indefinite_fluxes = np.count_nonzero(~np.isfinite(spectra), axis=0)

    number_fluxes = spectra.shape[0]
    keep_fluxes_mask = number_indefinite_fluxes < number_fluxes * drop

    spectra = spectra[:, keep_fluxes_mask]

    print(f"Shape after discard indefinite values: {spectra.shape}")

    # update wavelength grid
    wave = self.grid[keep_fluxes_mask]

    return [spectra, wave]
