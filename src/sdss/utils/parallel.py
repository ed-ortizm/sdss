"""Utility functionality for parallel computations"""

import numpy as np


def to_numpy_array(input_array, array_shape):
    """Create a numpy array backed by a shared memory Array."""

    share_array = np.ctypeslib.as_array(input_array)

    return share_array.reshape(array_shape)
