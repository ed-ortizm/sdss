"""Diferrent filters for spectra"""
from astropy.convolution import Gaussian1DKernel, convolve
import numpy as np

def filter_noise(spectrum, kernel_size):
    """Filter noise in spectra"""
    kernel = Gaussian1DKernel(kernel_size)

    extended_spectrum = np.hstack((spectrum[:10], spectrum, spectrum[-10:]))
    extended_spectrum = convolve(extended_spectrum, kernel)
    filtered_spectrum = extended_spectrum[10:-10]
    noise = spectrum - filtered_spectrum

    return noise, filtered_spectrum
