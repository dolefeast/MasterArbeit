import numpy as np
import scipy as sp

def remove_wiggles(x, y, p):
    """ Wiggle removing routine. Idea is to interpolate between the points in which the curvature goes to 0 to obtain the trend of the oscillation """
    x_ref = [0, 1]

    # Interpolate the given curve (x, y)
    # The trend of the curve should be where the curvature of the signal goes to 0
    fwiggle = sp.interpolate.UnivariateSpline(
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
    d2 = sp.interpolate.UnivariateSpline(x,
            derivs[2], 
            k=3,
            s=1.0
        )

    wzeros = d2.roots()
    wtrend = sp.interpolate.UnivariateSpline(
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
    return signal_filtered
