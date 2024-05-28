import numpy as np
def plot_from_0_to_1(array):
    """
    As the fields considered in my thesis are defined for z \in [0, 1],
    but the arrays don't know that, this is a way to reshape the arrays
    to carry that information
    """
    assert array.ndim == 1, "array must be 1-dimensional"
    z = list(range(len(array)))
    z = list(map(lambda x: x/len(array), z))

    return np.array(z), array
