import numpy as np
import scipy as sp


def modify_A0(z, rho):
    """
    Calculates the new vector field A0 corresponding to a family of solutions to the previous one
    Parameters:
        z [float]  the node array. An array of equally spaced z values from 0 to 1.
        rho [float] the rho at each node
    Returns:
        new_vector_field: callable, the electromagnetic field at z to be added to the base case
    """
    # Calculate this solving an ODE, not integrating

    # Need rho to be callable for solve_bvp
    try:
        iter(rho)  # Checking if rho is array
        rho_callable = sp.interpolate.CubicSpline(z, rho)
    except TypeError:
        try:
            rho(0.5)  # Check if rho is callable. z = 0.5 must be in the accepted values
            rho_callable = rho
            # if callable, do nothing
            pass
        except TypeError:  # rho was neither callable or array
            raise TypeError("rho should be either iterable or callable")

    def dAdz(z, A):
        return [A[1], -rho_callable(z)]

    def bc(Aa, Ab):
        return [
            Aa[1],
            Ab[0] - Aa[0],
        ]  # There must be no electric field on the boundary.

    factor = max(rho_callable(z))
    # rho_callable = lambda t, y: rho_callable(t)

    A0_guess = -factor * np.sin(2 * np.pi * z) / np.pi**2
    A0_gradient_guess = factor * np.cos(2 * np.pi * z) / np.pi

    modified_electric = sp.integrate.solve_ivp(
        lambda z, y: rho_callable(z), (0, 1), (0,), dense_output=True
    )

    electric_field = lambda t, y: -modified_electric.sol(t)

    modified_electric = sp.integrate.solve_ivp(
        # modified_electric.sol,
        electric_field,
        (0, 1),
        (0,),
        dense_output=True,
    )

    return modified_electric
