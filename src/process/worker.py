import numpy as np

###############################################################################
def to_numpy_array(shared_array, array_shape):
    """Create a numpy array backed by a shared memory Array."""
    share_array = np.ctypeslib.as_array(shared_array)
    return share_array.reshape(array_shape)


def init_worker(share_array, array_shape):

    global fluxes

    fluxes = to_numpy_array(share_array, array_shape)

    # global global_counter
    #
    # global_counter = counter
    #
    # # Update counter to show advance
    # with global_counter.get_lock():
    #     global_counter.value += 1
    #     counter_value = global_counter.value
    #     print(f"Get data from file N: {counter_value}")


def worker():
    """
    Robot computations for raw data
    """
    pass
