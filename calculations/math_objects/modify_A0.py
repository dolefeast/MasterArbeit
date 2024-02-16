import scipy as sp
import numpy as np


def modify_A0(self, charge_density):
    """
    Calculates the new vector field A0 corresponding to a family of solutions to the previous one
    Parameters:
        charge_density
    Returns:
        new_vector_field: callable, the electromagnetic field at z
    """
    # Calculate this solving an ODE, not integrating

    try:
        iter(charge_density)  # Checking if charge_density is array
        charge_density = sp.interpolate.CubicSpline(self.z, charge_density)
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

    print("self.n_points =", self.n_points)
    #print("np.shape(self.z) =", np.shape(self.n_points))
    print("np.shape(self.z) =", np.shape(self.z))

    modified_electric = sp.integrate.solve_bvp(
        dAdz, bc, self.z, np.zeros((self.n_points, 2)).T
    )

    print(np.shape(modified_electric.x))

    return modified_electric
