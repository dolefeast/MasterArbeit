import numpy as np


def iterable_output(func):
    def wrapper(*args):
        # Check which of the arguments are iterable
        result = []
        for i, arg in enumerate(args):
            try:
                iter(arg)  # arg is the array, the others are fixed

                for arg_value in arg:
                    all_scalar_arg = list(args).copy()
                    all_scalar_arg[i] = arg_value
                    result.append(func(*all_scalar_arg))
                return np.array(result)
            except TypeError:
                return func(*args)  # All inputs are scalar

    return wrapper
