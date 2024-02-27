import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scripts.filtering import moving_average
from scripts.plotting import plot_from_0_to_1

def check_linear(x):
    return all(np.diff(x)==np.diff(x)[0])

def remove_noisy_neighbourhood(
        x, 
        y,
        points=[float],
        size=float):
    """
    Given curve (x, y) with problematic points=(x1, x2, ...), take their neighbourhood with size=size away and interpolate around it, thus smoothing the curve. 
    """
    x_list = list(x)
    y_list = list(y)
    x_array = np.array(x)
    idx = []
    # The points to take out of the array
    for p in points:
        idx += list(np.where(abs(x_array - p)<=size/2))
    idx = np.reshape(
                np.array(idx),
                -1 # 1-d array
            )
    x_list = [x for i, x in enumerate(x_list) if i not in idx]
    y_list = [y for i, y in enumerate(y_list) if i not in idx]

    return np.array(x_list), np


def periodic_signal(x, y, padding_size=None):
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


def remove_wiggles(x, y):
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

    return wtrend

def main():
    unfiltered_charge_density = np.genfromtxt('./saved_solutions/lambda_1_mass_1.csv', delimiter=",")

    x, y = unfiltered_charge_density
    x, y= periodic_signal(x, y)

    no_wiggle_y = remove_wiggles(x, y)
    plt.plot(x, y, alpha= 0.4)
    plt.plot(x, no_wiggle_y(x))

    plt.show()

if __name__ == "__main__":
#    main()
    x = np.linspace(-1, 2, 500)
    #plt.plot(x)
    x,y = remove_noisy_neighbourhood(x, x, [0,1], 0.3)

    plt.plot(x, x, 'o')
    plt.show()
