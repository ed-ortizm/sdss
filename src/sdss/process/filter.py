"""Diferrent filters for spectra"""
from astropy.convolution import Gaussian1DKernel, convolve
import numpy as np


def filter_noise(spectrum: np.array, kernel_size: int) -> tuple:
    """
    Filter noise on a spectrum with a gaussian kernel

    INPUT

    spectrum: data with noise
    kernel_size: number of elements in gaussian kernel

    OUTPUT

    filtered_spectrum, noise:
        filtered_spectrum: spectrum with noise removed
        noise: spectrum's noise
    """
    kernel = Gaussian1DKernel(kernel_size)

    filtered_spectrum = convolve(spectrum, kernel, boundary="extend")
    noise = spectrum - filtered_spectrum

    return filtered_spectrum, noise
