import scipy as sp


def normalize(field: callable):
    """
    Returns the normalized field
    Parameter:
       field: callable
            The 'quantum' field to be integrated over
    Returns:
        normalized_field: callable
            The normalized 'quantum' field
    """
    norm_squared = sp.integrate.quad(lambda z: field(z) * sp.conjugate(field(z)), 0, 1)

    return lambda z: field(z) / np.sqrt(norm_squared)
