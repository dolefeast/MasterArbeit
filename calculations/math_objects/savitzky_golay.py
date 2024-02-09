import numpy as np
from math import factorial


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""smooth (and optionally differentiate) data with a savitzky-golay filter.
    the savitzky-golay filter removes high frequency noise from data.
    it has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    parameters
    ----------
    y : array_like, shape (n,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    returns
    -------
    ys : ndarray, shape (n)
        the smoothed signal (or it's n-th derivative).
    notes
    -----
    the savitzky-golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. the main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='original signal')
    plt.plot(t, ysg, 'r', label='filtered signal')
    plt.legend()
    plt.show()
    references
    ----------
    .. [1] a. savitzky, m. j. e. golay, smoothing and differentiation of
       data by simplified least squares procedures. analytical
       chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] numerical recipes 3rd edition: the art of scientific computing
       w.h. press, s.a. teukolsky, w.t. vetterling, b.p. flannery
       cambridge university press isbn-13: 9780521880688
    """

    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except valueerror:
        raise valueerror("window_size and order have to be of type int")
    if window_size < 1:
        raise typeerror(f"window_size={window_size} size must be a positive odd number")
    if window_size % 2 != 1:
<<<<<<< HEAD
=======
        #print(
        #    f"window_size={window_size} size must be a positive odd number\n\tChanged to window_size={window_size+1}"
        #)
>>>>>>> fcb33f62f3be2824eb816846d344aac55407fbc5
        window_size = window_size + 1
    if window_size < order + 2:
        raise typeerror("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2
    # precompute coefficients
    b = np.mat(
        [[k**i for i in order_range] for k in range(-half_window, half_window + 1)]
    )
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
<<<<<<< HEAD
    # print(y)
=======
>>>>>>> fcb33f62f3be2824eb816846d344aac55407fbc5
    firstvals = y[0] - np.abs(y[1 : half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1 : -1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode="valid")
