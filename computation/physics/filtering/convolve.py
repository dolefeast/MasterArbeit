import numpy as np

def convolve(
        self,
    x: [float],
    y: [float],
):
    """
    Convolves the signal with a very specific function to smooth it out.
    Parameters:
        x: [float], the x array #Not used, but it is needed for consistency
        y: [float], the signal
        delta_N: int, the size of the convolution array. MUST be equal to
            N_points/(N_mode_cutoff+1)
    Returns:
        x, y_filtered: ([float], [float]), (the original x-array, the filtered signal)
    ---------------
    The very specific array is handpicked to cancel out the oscillations
    that modulate the curve. It is of the shape [3, ..., 0.4, ..., 6, ..., 0.4, ..., 3]
    where its length is delta_N=self.n_points/(self.N_mode_cutoff + 1). delta_N is the size
    of a whole oscillation of the noise to be filtered out.
    """

    delta_N = self.n_points // ( self.N_mode_cutoff  + 1) + 1
    # +1 since otherwise the last mesh point is not reached 
    array = [0] * delta_N

    surround_peaks = 0.00
    peaks = 0.4
    middle = 6
    edges = 3

    # On the borders
    array[0] = edges
    array[-1] = edges

    # Peak of the sine
    array[delta_N // 4 + 1] = surround_peaks
    array[delta_N // 4 - 1] = surround_peaks

    # Trough of the sine
    array[3 * (delta_N // 4) + 1] = surround_peaks
    array[3 * (delta_N // 4) - 1] = surround_peaks

    array[delta_N // 4] = peaks
    array[3 * (delta_N // 4)] = peaks

    # In the middle of the curve
    array[delta_N // 2] = middle

    # To normalize
    array = np.array(array) / sum(array)

    y_filtered = np.convolve(y, array, mode="same")

    return x, y_filtered

def convolve_twice(
        self,
    x: [float],
    y: [float],
    ):
    """
    Convolves the signal twice with a very specific function to smooth it out.
    Parameters:
        x: [float], the x array #Not used, but it is needed for consistency
        y: [float], the signal
        delta_N: int, the size of the convolution array. MUST be equal to
            N_points/(N_mode_cutoff+1)
    Returns:
        y_filtered: [float], the filtered signal
    ---------------
    The very specific array is handpicked to cancel out the oscillations
    that modulate the curve. It is of the shape [3, ..., 0.4, ..., 6, ..., 0.4, ..., 3]
    where its length is delta_N=self.n_points/(self.N_mode_cutoff + 1). delta_N is the size
    of a whole oscillation of the noise to be filtered out.
    """
    x, y_filtered = self.convolve(
            x, 
            y,
            )

    x, y_filtered_twice = self.convolve(
            x, 
            y_filtered,
            )

    return x, y_filtered_twice
