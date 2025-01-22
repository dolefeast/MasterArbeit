import numpy as np
from scipy.interpolate import CubicSpline

def extend_signal(
        self,
        x,
        y,
        padding_size=None
        ):
    assert len(x) == len(y), "x and y don't have the same shape"

    # if padding_size is not specified,
    # then the padding is y itself
    if padding_size is None:
        padding_size = len(y)
    # wrap gives periodization of the function
    y_padded = np.pad(y, padding_size, mode="wrap")

    # need to extend x, too
    dx = x[1] - x[0]

    x_padded = np.linspace(-dx * padding_size, 1 + dx * padding_size, len(y_padded))

    return x_padded, y_padded

def remove_neighbourhood(
        self,
    x,
    y,
    points: (float)=(0,1),
    neighbourhood_size: float=None,
    force_zero: bool=True,
):
    """
    Given curve (x, y) with problematic neighbourhoods around points=(x1, x2, ...), take their neighbourhood with neighbourhood_size=neighbourhood_size away and interpolate around it, thus smoothing the curve.
    """

    if neighbourhood_size is None:
        # Take out the points corresponding to the convolution array
        neighbourhood_size = (
        1 / ( self.max_N + 1 ) 
        + 1 / self.n_points # Without this, the algorithm takes away
                            # an "open" neighbourhood (can't be open since it's discrete), 
                            # and adding it converts it into a closed one.
                            # Corresponds to adding a dz.
            )

    x_list = list(x)
    y_list = list(y)
    
    x_array = np.array(x)
    
    # A python array to use its addition properties
    idx = []

    # The points to take out of the array
    for p in points:
        
        window = np.where(abs(x_array - p) <= neighbourhood_size)[0].tolist()
        idx += window

    idx = np.reshape(
        np.array(idx),
        -1,  # 1-d array
    )
    x_list = [x for i, x in enumerate(x_list) if i not in idx]
    y_list = [y for i, y in enumerate(y_list) if i not in idx]

    if force_zero:  # Force the removed values to go through 0
        for p in points:
            for x_index, x_value in enumerate(x_list):
                if x_value > p:
                    x_list = x_list[:x_index] + [p] + x_list[x_index:]
                    y_list = y_list[:x_index] + [0] + y_list[x_index:]
                    break

    return np.array(x_list), np.array(y_list)

def remove_and_interpolate(
        self,
    x: [float],
    y: [float],
    points=(0, 1),
    neighbourhood_size:float=None,
    force_zero: bool=True,
):
    x_removed, y_removed = self.remove_neighbourhood(
        x,
        y,
        points=points,
        neighbourhood_size=neighbourhood_size,
        force_zero=True
    )

    interpolated_curve = CubicSpline(
            x_removed,
            y_removed,
            )

    return x, interpolated_curve(x)  # So that both are arrays

def return_to_0_1(
        self,
        x, 
        y):

    idx = np.where(
        np.logical_and(
            x >= 0,
            x <= 1,
        )
    )
    return x[idx], y[idx]

def filter_rho( self,rho):
    """
    Parameters:
        rho: a noisy signal (with a very particular kind of noise)
    Returns:
        rho_filtered: the filtered signal.
    """
    # This exploits the fact that we know exactly what the 
    # noise looks like. It can be convoluted away with
    # a very specific function. Not particularly
    # relevant.

    # As long as n_points is exactly 8 * (max_N + 1),
    # the filtering should work fine
    # As with A0, n_points, max_N are necessarily
    # defined somewhere before the function is called.
    delta_N = self.n_points // ( self.max_N + 1) + 1
    convoluting_array = [0] * delta_N

    # Some fine tuned parameters
    surround_peaks = 0.00
    peaks = 0.4
    middle = 6
    edges = 3

    # On the borders
    convoluting_array[0] = edges
    convoluting_array[-1] = edges

    # Peak of the sine
    convoluting_array[delta_N // 4 + 1] = surround_peaks
    convoluting_array[delta_N // 4 - 1] = surround_peaks

    # Trough of the sine
    convoluting_array[3 * (delta_N // 4) + 1] = surround_peaks
    convoluting_array[3 * (delta_N // 4) - 1] = surround_peaks

    convoluting_array[delta_N // 4] = peaks
    convoluting_array[3 * (delta_N // 4)] = peaks

    # In the middle of the curve
    convoluting_array[delta_N // 2] = middle

    # Normalizing it so that we don't get extra factors
    convoluting_array = np.array(convoluting_array) / sum(convoluting_array)

    rho_filtered = np.convolve(rho, convoluting_array, mode="same")
 
    return rho_filtered # This will stield yield some noise as 
                      # this is not the full smoothing algorithm.
def extend_and_filter(
        self,
    x,
    y,
    neighbourhood_size=None,
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
    x_extended, y_extended = self.extend_signal(x, y, padding_size=padding_size)

    # Second filter it
    y_extended_filtered = self.filter_rho(y_extended)
    y_extended_filtered = self.filter_rho(y_extended_filtered)

    # Third and fourth remove the boundaries and interpolate
    x_extended_removed, y_extended_filtered_removed = self.remove_and_interpolate(
        x_extended,
        y_extended_filtered,
        points=points,
        neighbourhood_size=neighbourhood_size,
        force_zero=force_zero,
    )

    # Fifth return the signal only in the interval 0 to 1
    x_0_to_1, y_0_to_1 = self.return_to_0_1(x_extended_removed, y_extended_filtered_removed)

    return x_0_to_1, y_0_to_1
