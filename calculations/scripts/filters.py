import numpy as np
import scipy as sp
import scipy.interpolate as interpolate

def convolve(
        x: [float],
        y: [float],
        delta_N=9,
    ):
    """ 
    Convolves the signal with a very specific function to smooth it out.
    Parameters:
        x: [float], the x array #Not used, but it is needed for consistency
        y: [float], the signal  
        delta_N: int, the size of the convolution array. MUST be equal to
            N_points/(N_mode_cutoff+1)
    Returns:
        y_filtered: [float], the filtered signal
    ---------------
    The very specific array is handpicked to cancel out the oscillations 
    that modulate the curve. It is of the shape [3, 0, 0.4, 0, 6, 0, 0.4, 0, 3] 
    and if delta_N doesn't fit exactly, the filtering does not work
    If delta_N changes, it just padds with zeroes.
    """

    delta_N = 9
    array = [0] * delta_N

    surround_peaks = 0.00
    peaks = 0.4
    middle = 6
    edges = 3

    #On the borders
    array[0] = edges
    array[-1] = edges

    # Peak of the sine
    array[delta_N//4+1] = surround_peaks
    array[delta_N//4-1] = surround_peaks

    # Trough of the sine
    array[3 * ( delta_N//4 )+1] = surround_peaks
    array[3 * ( delta_N//4 )-1] = surround_peaks

    array[delta_N//4] = peaks
    array[3 * ( delta_N//4 )] = peaks

    # In the middle of the curve
    array[delta_N//2] = middle

    #To normalize
    array = np.array(array) / sum(array)

    print(y, array)

    y_filtered = np.convolve(y, array, mode='same')

    return y_filtered

def convolve_twice(x, y, delta_N):
    """
    Applies convolve twice
    """
    y_filtered = convolve(x, y, delta_N)
    return convolve(x, y_filtered, delta_N)

def remove_wiggles(x, y, p):
    """ Wiggle removing routine. Idea is to interpolate between the points in which the curvature goes to 0 to obtain the trend of the oscillation """
    x_ref = [0, 1]

    # Interpolate the given curve (x, y)
    # The trend of the curve should be where the curvature of the signal goes to 0
    fwiggle = interpolate.UnivariateSpline(
            x,
            y,
            k=3,
            s=0
            )
    # fwiggle is a cubic spline. 
    # derivs is an array which contains the 3 derivatives (0, 1 and 2 order) at each point of x
    derivs = np.array(
        [fwiggle.derivatives(_k) for _k in x]
            ).T
    # Get and interpolate the second derivative
    d2 = interpolate.UnivariateSpline(x,
            derivs[2], 
            k=3,
            s=1.0
        )

    wzeros = d2.roots()
    wtrend = interpolate.UnivariateSpline(
            wzeros,
            fwiggle(wzeros),
            k=3,
            s=0
            )

    return wtrend(x)

def remove_wiggles_twice(x, y, p):
    y = remove_wiggles(x, y, p)
    return remove_wiggles(x, y, p)

def remove_wiggles_thrice(x, y, p):
    y = remove_wiggles_twice(x, y, p)
    return remove_wiggles(x, y, p)

def moving_average(x, signal, w):
    signal_filtered = np.convolve(signal, np.ones(w), 'full') / w
    X = np.linspace(
               x[0],
               x[-1],
               len(signal_filtered
                )
               )

    signal_interpolated = interpolate.UnivariateSpline(
        X,
            signal_filtered,
            k=3,
            s=0
            )
    signal_x_shape = signal_interpolated(x)
    return signal_x_shape

def moving_average_twice(x, signal, w):
    signal_filtered = np.convolve(signal, np.ones(w), 'full') / w
    signal_interpolated = interpolate.UnivariateSpline(
            x, 
            signal_filtered,
            k=3,
            s=0
            )
    signal_x_shape = signal_interpolated(x)
    return moving_average(x, signal_filtered, w)

