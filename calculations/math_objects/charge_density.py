from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.savitzky_golay import savitzky_golay
from math_objects.normalize import normalize
from math_objects.get_eigenvalues import get_eigenvalues

import scipy as sp
import numpy as np
import mpmath

def normalized_eigenstate(self,
    eigenstate,
    integration_limit=200):
    """
    Calculates the charge density associated to certain eigenstate
    Parameters:
        eigenstate: scipy.integrate.solve_bvp solution object
        integration_limit: float = 200. The integration limit for scipy.integrate.quad in the normalization condition
    Returns:
        charge_density: callable, the charge density corresponding to the input eigenstate
    """

    eigenvalue = eigenstate.p[0]
    charge_density_without_normalization = lambda z: (eigenvalue - self.e * self.A0_field(z)) * np.abs(eigenstate.sol(z)[0])**2

    # calculating its norm
    charge_density_norm_squared = mpmath.quad(
        charge_density_without_normalization, [0, 1]
    )

    charge_density_normalized = lambda z:  np.sign(eigenvalue)*(
        charge_density_without_normalization(z) 
        / charge_density_norm_squared
        )

    return charge_density_normalized

def calculate_total_charge_density(self,
        eigenstates,
        renormalization=False,
        smoothing:str='savgol'):

    """
    Calculates total charge density as a callable from an array of eigenstates
    Parameters:
        eigenstates: [scipy.integrate.solve_bvp] n-array of solution objects
        renormalization: False. Whether renormalization is applied (explained below)
        filtering: str='savgol'. Can be False. Whether filtering is applied to each eigenstate, 
                                 and which filtering is applied.
    Returns:
        total_charge_density: callable, the charge density corresponding to the input eigenstates
    """
    
    def total_charge_density(z, renormalization=renormalization):
        n_index_array = np.sqrt(get_eigenvalues(eigenstates)**2 - self.m**2)/np.pi
        if renormalization:
            renormalization_closed_form = (
                    self.e
                    * self.lambda_value 
                    * z * (1-z)
                    / np.tan(np.pi*z))
        else:
            renormalization_closed_form = 0


        n_index_array = n_index_array

        adding_sign = 1 
        true_eigenstate_array = eigenstates
        total_charge_density = 0
        for eigenstate, n in zip(eigenstates, n_index_array):
            if renormalization:
                renormalization_term = (
                              adding_sign
                            * self.e 
                            * self.lambda_value 
                            * z * (1-z)
                            * np.sin(np.pi * n * z)
                            * np.cos(np.pi * n * z)
                            )
            else:
                renormalization_term = 0

            total_charge_density += (
                    0.5*self.normalized_eigenstate(eigenstate)(z) 
                    + renormalization_term
                                )
            potential_term = self.e**2 / np.pi * self.A0_field(z)

        total_charge_density -= adding_sign * renormalization_closed_form 
        total_charge_density *= self.e

        if smoothing:
            total_charge_density = savitzky_golay(
                    total_charge_density,
                    int(90/1000 * self.n_points),
                    3
                    )

        return (
                total_charge_density
                 + potential_term 
                )

    return total_charge_density


    
