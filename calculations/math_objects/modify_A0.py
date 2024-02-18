import scipy as sp
import numpy as np

def modify_A0(z, charge_density):
    """
    Calculates the new vector field A0 corresponding to a family of solutions to the previous one
    Parameters:
        z [float]  the node array. An array of equally spaced z values from 0 to 1.
        charge_density [float] the charge_density at each node
    Returns:
        new_vector_field: callable, the electromagnetic field at z to be added to the base case
    """
    # Calculate this solving an ODE, not integrating

    # Need charge_density to be callable for solve_bvp
    try: 
        iter(charge_density)  # Checking if charge_density is array
        charge_density = sp.interpolate.CubicSpline(z, charge_density)
    except TypeError:
        try:
            charge_density(
                0.5
            )  # Check if charge_density is callable. z = 0.5 must be in the accepted values
        # if callable, do nothing
            pass
        except TypeError:
            raise TypeError("charge_density should be either iterable or callable")

    def dAdz(z, A):
        return [A[1], -charge_density(z)]

    def bc(Aa, Ab):
        return [Aa[0], Ab[0]]  # There must be no potential on the boundary.


    modified_electric = sp.integrate.solve_bvp(
        dAdz, bc, z, np.zeros_like((z, z))
    )

    return modified_electric
