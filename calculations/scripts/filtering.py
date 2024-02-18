import numpy as np

def moving_average(signal, w):
    return np.convolve(signal, np.ones(w), 'full') / w

def double_filtering(
        signal: [float],
        window: [int],
        ):
    first_filter = moving_average(signal, window//4+1)
    second_filter = moving_average(first_filter, window)
    return second_filter

#TODO
def recursive_filtering(
        signal: [float],
        filter_method: callable,
        filter_parameters: tuple, 
        n_recursion: int
        ):
    """ Recursively applies filter
    Parameters:
        signal: list. 
        filter_method: callable. Call signal: filter_method(signal, *filter_parameters) The filter method to be applied
        filter_parameters. the parameters to be passed to the filtering method
        n_recursion: int=3: How many times the filtering is applied
    returns:
        filtered_signal: [float]
    """
    # When n_recursion reaches 0
    print('The filter_parameters are:', filter_parameters)
    if not n_recursion:
        return filter_method(signal, *filter_parameters)
    return recursive_filtering(signal, filter_method, *filter_parameters[:-1], filter_parameters-1)
