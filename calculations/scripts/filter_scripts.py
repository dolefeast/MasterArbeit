import numpy as np
import scipy as sp

def extend_signal(x, y, padding_size=None):
    assert len(x) == len(y)
    
    # if padding_size is not specified,
    # then the padding is y itself
    if padding_size is None:
        padding_size = len(y)
    # wrap gives periodization of the function
    y_padded = np.pad(y, padding_size, mode='wrap') 

    # need to periodize x
    dx = x[1] - x[0]
    
    x_padded = np.linspace(
            -dx*padding_size,
            1+dx*padding_size,
            len(y_padded)
            )

    return x_padded, y_padded

def remove_neighbourhood(
        x, 
        y,
        points:[float],
        neighbourhood_size:float,
        force_zero: bool=True):
    """
    Given curve (x, y) with problematic points=(x1, x2, ...), take their neighbourhood with neighbourhood_size=neighbourhood_size away and interpolate around it, thus smoothing the curve. 
    """
    x_list = list(x)
    y_list = list(y)
    x_array = np.array(x)
    idx = []
    # The points to take out of the array
    for p in points:
        window = np.where(abs(x_array - p)<=neighbourhood_size)[0].tolist()
        idx += window

    idx = np.reshape(
                np.array(idx),
                -1 # 1-d array
            )
    x_list = [x for i, x in enumerate(x_list) if i not in idx]
    y_list = [y for i, y in enumerate(y_list) if i not in idx]

    if force_zero: # Force the removed values to go through 0
        for p in points:
            for x_index, x_value in enumerate(x_list):
                if x_value > p:
                    x_list = x_list[:x_index] + [p] + x_list[x_index:]
                    y_list = y_list[:x_index] + [0] + y_list[x_index:]
                    break

    return np.array(x_list), np.array(y_list)

def remove_and_interpolate(
        x: [float],
        y: [float],
        points=(0,1),
        neighbourhood_size=float,
        force_zero:bool=True,
        ):
    x_removed, y_removed = remove_neighbourhood(x, 
            y, 
            points=points, 
            neighbourhood_size=neighbourhood_size,
            force_zero=True)

    interpolated_curve = sp.interpolate.UnivariateSpline(
            x_removed, 
            y_removed,
            k=3,
            s=0
            )
    return x, interpolated_curve(x) # So both are arrays

def return_to_0_1(x, y):
    idx = np.where(np.logical_and
            (
                x>=0,
                x<=1,
                )
            )
    return x[idx], y[idx]


def extend_and_filter(
    x,
    y,
    filter_method,
    filter_parameters,	
    neighbourhood_size,
    padding_size=None,
    points=(0, 1),
    force_zero=True,
    ):
    """
    Parameters:
    x: [float], the x-array of the signal. increasing and goes from 0 to 1
    y: [float], the y-array of the signal
    filter_method: callable, the method to be used for filtering
    neighbourhood_size: float, the 'diameter' of the neighbourhood to be removed around points
    padding_size=None, the size of the padding i.e. the extension of the signal. if none, padding_size = len(y) and therefore the signal is copied both ways
    points=(0, 1), the points around which the neighbourhood is to be removed
    force_zero:bool=True, whether after removing a certain neighbourhood, we force y(points) to go through 0
    filter_parameters=(, ),	 filter parameters that the filtering script may need

    1. Extends signal by padding_size
    2. Filters the extended signal
    3. Removes the points on a neighbourhood of the boundaries (0, 1)
     3a. If force_values=True, force the signal to go through the forced values (initially 0)
    4. Interpolates over the void area
    5. Returns the signal only in the interval (0, 1)

    """
    # First extend the signal
    x_extended, y_extended = extend_signal(x, y, padding_size=padding_size)

    # Second filter it
    y_extended_filtered = filter_method(
            x_extended,
            y_extended,
            *filter_parameters
            )

    # Third and fourth remove the boundaries and interpolate
    x_extended_removed, y_extended_filtered_removed = remove_and_interpolate(
            x_extended,
            y_extended_filtered,
            points=points,
            neighbourhood_size=neighbourhood_size,
            force_zero=force_zero,
            )

    # Fifth return the signal only in the interval 0 to 1
    x_0_to_1, y_0_to_1 = return_to_0_1(x_extended_removed, y_extended_filtered_removed)

    return x_0_to_1, y_0_to_1

def double_filtering(
        signal: [float],
        window: [int],
        ):
    first_filter = moving_average(signal, window//4+1)
    second_filter = moving_average(first_filter, window)
    return second_filter

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
