import numpy as np
###############################################################################
def init_worker(counter):
    """
    Initialize worker for raw data process
    PARAMETERS
        counter:
    """
    global global_counter

    # Update counter to show advance
    with global_counter.get_lock():
        global_counter.value += 1
        counter_value = global_counter.value

    print(f"Get data from file N: {counter_value}", end="\r")


def worker():
    """
    Robot computations for raw data
    """
    pass
