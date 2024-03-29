import numpy as np
import scipy as sp


def extend_signal(x, y, padding_size=None):
    assert len(x) == len(y)

    # if padding_size is not specified,
    # then the padding is y itself
    if padding_size is None:
        padding_size = len(y)
    # wrap gives periodization of the function
    y_padded = np.pad(y, padding_size, mode="wrap")

    # need to periodize x
    dx = x[1] - x[0]

    x_padded = np.linspace(-dx * padding_size, 1 + dx * padding_size, len(y_padded))

    return x_padded, y_padded


def remove_wiggles(x, y):
    """Wiggle removing routine. Idea is to interpolate between the points in which the curvature goes to 0 to obtain the trend of the oscillation"""
    x_ref = [0, 1]

    # Interpolate the given curve (x, y)
    # The trend of the curve should be where the curvature of the signal goes to 0
    fwiggle = sp.interpolate.UnivariateSpline(x, y, k=3, s=0)
    # fwiggle is a cubic spline.
    # derivs is an array which contains the 3 derivatives (0, 1 and 2 order) at each point of x
    derivs = np.array([fwiggle.derivatives(_k) for _k in x]).T
    # Get and interpolate the second derivative
    d2 = sp.interpolate.UnivariateSpline(x, derivs[2], k=3, s=1.0)

    wzeros = d2.roots()
    wtrend = sp.interpolate.UnivariateSpline(wzeros, fwiggle(wzeros), k=3, s=0)

    return wtrend(x)


def remove_neighbourhood(x, y, points: [float], size: float, force_zero: bool = True):
    """
    Given curve (x, y) with problematic points=(x1, x2, ...), take their neighbourhood with size=size away and interpolate around it, thus smoothing the curve.
    """
    x_list = list(x)
    y_list = list(y)
    x_array = np.array(x)
    idx = []
    # The points to take out of the array
    for p in points:
        window = np.where(abs(x_array - p) <= size / 2)[0].tolist()
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
                    x[x_index:x_index] = p
                    y[x_index:x_index] = 0

    return np.array(x_list), np.array(y_list)


def remove_and_interpolate(
    x: [float],
    y: [float],
    points=(0, 1),
    size=float,
):
    x_holes, y_holes = remove_neighbourhood(x, y, points=points, size=size)
    interpolated_curve = sp.interpolate.UnivariateSpline(x_holes, y_holes, k=3, s=0)
    return x, interpolated_curve(x)  # So both are arrays


def return_to_0_1(x, y):
    idx = np.where(
        np.logical_and(
            x >= 0,
            x <= 1,
        )
    )
    return x[idx], y[idx]


def bao_filtering(
    x: [float], y: [float], points=(0, 1), size: float = 0.01, padding_size=None
):
    """
    Filters the signal based on BAO signal filtering.
    Parameters
        x: [float], # Array from 0 to 1
        y: [float], # Signal to be filtered
        points=(0,1), # Problematic points to be taken out and interpolated through
        size=float, # Size of the window to be taken out.
        padding_size=None # How much should the y array be extended. If None,
                          # y gets extended by y in both directions
    Returns:
        x: [float]  #Array from 0 to 1
        y: [float]  #Filtered signal

    The routine is as follows:
        1. extend_signal: Extend the signal periodically to the left of 0 and right of 1.
            Needed because otherwise problems with the derivative of the filtered signal at the boundaries.
        2. remove_wiggles: Interpolate curve, calculate its second derivative, calculate the zeros. Interpolate between the zeros of the second derivative.
        3. remove_and_interpolate: This gives noise at 0 and 1 boundaries. Remove these points and interpolate.
        4. return_to_0_1: filtered signal has domain (0-padding_size, 1+padding_size). Return to domain (0, 1).
    """
    x_extended, y_extended = extend_signal(x, y, padding_size=None)
    no_wiggle_y = remove_wiggles(x_extended, y_extended)

    x_no_boundaries, y_no_boundaries = remove_and_interpolate(
        x_extended, no_wiggle_y, points=points, size=size
    )

    return return_to_0_1(x_no_boundaries, y_no_boundaries)
