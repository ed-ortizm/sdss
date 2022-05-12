"""Module to measure execution time of functions"""
from functools import wraps
import time


def timeit(function):
    """Measure execution time of function"""

    @wraps(function)
    def timeit_wrapper(*args, **kwargs):

        start_time = time.perf_counter()
        result = function(*args, **kwargs)
        finish_time = time.perf_counter()

        print(f"Function {function.__name__} took {finish_time-start_time}")

        return result

    return timeit_wrapper
