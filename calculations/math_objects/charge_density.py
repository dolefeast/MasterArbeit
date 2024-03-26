from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.get_eigenvalues import get_eigenvalues

import scipy as sp
import numpy as np
import mpmath

def normalized_eigenstate(self,
    eigenstate,
    integration_limit=300):
    """
    Calculates the charge density rho associated to certain eigenstate
    Parameters:
        eigenstate: scipy.integrate.solve_bvp solution object
        integration_limit: float = 200. The integration limit for scipy.integrate.quad in the normalization condition
    Returns:
        rho: callable, the charge density corresponding to the input eigenstate
    """

    eigenvalue = eigenstate.p[0]
    rho_without_normalization = lambda z: (eigenvalue - self.e * self.A0_field(z)) * np.abs(eigenstate.sol(z)[0])**2

    # calculating its norm
    rho_norm_squared = mpmath.quad(
        rho_without_normalization, [0, 1]
    )

    rho_normalized = lambda z:  np.sign(eigenvalue)*(
        rho_without_normalization(z) 
        / rho_norm_squared
        )

    return rho_normalized

def calculate_total_rho(self,
        renormalization=False,
        ):

    """
    Calculates total charge density as a callable from an array of eigenstates
    Parameters:
        eigenstates: [scipy.integrate.solve_bvp] n-array of solution objects
        renormalization: False. Whether renormalization is applied (explained below)
        filtering: str='savgol'. Can be False. Whether filtering is applied to each eigenstate, 
                                 and which filtering is applied.
    Returns:
        total_rho: callable, the charge density corresponding to the input eigenstates
    """
    
    def total_rho(z, renormalization=renormalization):
        N = len(self.eigenstate_array)
        n_index_array = np.arange(-N//2, N//2+1)
        if renormalization:
            renormalization_closed_form = (
                    self.e
                    * self.lambda_value 
                    * z * (1-z)
                    / np.tan(np.pi*z))
        else:
            renormalization_closed_form = 0


        adding_sign = 1 
        true_eigenstate_array = self.eigenstate_array
        total_rho = 0
        for n, rho_n in zip(n_index_array, self.rho_n_array): 
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

            total_rho += (
                    0.5*rho_n
                    + renormalization_term
                                )
            potential_term = self.e**2 / np.pi * self.A0_field(z)

        total_rho -= adding_sign * renormalization_closed_form 
        total_rho *= self.e

        return (
                total_rho
                 + potential_term 
                )

    return total_rho


    
