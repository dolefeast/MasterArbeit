import numpy as np
from itertools import count
def value_generator(starting_value, ratio, method='linear'):
    """
    An element generator to yield values in different behaviours.
    Parameters:
        starting_value: float. The value at which the generator will start
        ratio: float. some other parameter to measure the speed at which the steps are taken
        method: string or f(x, ratio). if a string, should be either 'linear', 'quadratic', 'logarithmic', 'exponential'
    """
    methods = {'linear': lambda x, ratio: ratio*x, 
            'quadratic': lambda x, ratio: ratio * x ** 2,
            'logarithmic': lambda x, ratio: ratio * np.log(1+x) , #So that it starts in 0, not 1
            'exponential': lambda x, ratio: ratio * np.exp(x)}

    if method not in methods:
        if isinstance(method, str):
            raise ValueError(f"The specified method {method} is not in {methods.keys()}.")
        elif callable(method):
            try:
                method(0, 0)
                f = method
            except TypeError:
                raise TypeError("The supplied function is not of the form f(starting_value, ratio)")
    else:
        f = methods[method]

    for i in count():
        yield f(i+starting_value, ratio)

