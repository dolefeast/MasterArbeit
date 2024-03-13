import numpy as np
import scipy as sp

import matplotlib.pyplot as plt


from mpmath import odefun

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
        charge_density_callable = sp.interpolate.CubicSpline(z, charge_density)
    except TypeError:
        try:
            charge_density(
                0.5
            )  # Check if charge_density is callable. z = 0.5 must be in the accepted values
            charge_density_callable = charge_density
        # if callable, do nothing
            pass
        except TypeError: # charge_density was neither callable or array
            raise TypeError("charge_density should be either iterable or callable")

    def dAdz(z, A):
        return [A[1], -charge_density_callable(z)]

    def bc(Aa, Ab):
        return [Aa[1], Ab[0]-Aa[0]]  # There must be no electric field on the boundary.


    factor = max(charge_density_callable(z))
    # charge_density_callable = lambda t, y: charge_density_callable(t)

    A0_guess = -factor*np.sin(2*np.pi*z)/np.pi**2
    A0_gradient_guess = factor*np.cos(2*np.pi*z)/np.pi

    modified_electric = sp.integrate.solve_ivp(
            lambda z, y: charge_density_callable(z),
            (0, 1),
            (0, ),
            dense_output=True
            )
    
    electric_field = lambda t,y: -modified_electric.sol(t)

    modified_electric = sp.integrate.solve_ivp(
            # modified_electric.sol,
            electric_field,
            (0, 1),
            (0, ),
            dense_output=True
            )

    return modified_electric
