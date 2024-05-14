def total_filtering_dirichlet(
        self,
        z, 
        signal,
        ):
    """
    Parameters:
        z, the independent variable. Should be equal to self.z. 
        signal, the function of z to be filtered,
    Returns:
        z, of the original shape
        signal_filtered, the filtered signal (assuming dirichlet)
    ----------------------------
    Given a certain noisy signal (whose noise function has a period = 1 / (self.N_mode_cutoff + 1))
    filter the signal following the algorithm:
        1. Extend the signal periodically both to the left and to the right
            to preserve the derivative at the boundaries
        2. Filter it using a hand picked convolution array.
        3. Since the signal won't behave nicely at the boundaries,
            take them away, force the function to go through 0 at z=0, 1.
            Then interpolate
        4. Return to the original shape of the signal i.e. with z = [0, 1]
    """

    # First extend
    z, signal = self.extend_signal(z, signal)

    # Then filter
    z, signal = self.convolve_twice(
            z,
            signal
            )

    # Take out 0 1 noisy neighbousignalods
    z, signal = self.remove_and_interpolate(
            z,
            signal
            )

    # Return to z \in [0, 1]
    z, signal = self.return_to_0_1(
            z,
            signal
            )

    return z, signal
