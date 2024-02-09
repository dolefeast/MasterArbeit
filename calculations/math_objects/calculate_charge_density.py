import numpy as np
import scipy as sp

def calculate_charge_density(self, solution):
    """
    Calculates the charge density associated to certain eigenstate
    Parameters:
        eigenstates: scipy.integrate.solve_bvp solution object
    Returns:
        charge_density: callable, the charge density corresponding to the input eigenstate
    """

    eigenvalue = solution.p[0]
    charge_density_no_normalization = lambda z: (
        eigenvalue - self.phi.charge * self.A0(z)
    ) * np.real(solution.sol(z)[0] ** 2)
    charge_density_norm_squared = sp.integrate.quad(
        charge_density_no_normalization, 0, 1, limit=200
    )[0]
    charge_density_normalized = (
        lambda z: np.sign(eigenvalue)
        * charge_density_no_normalization(z)
        / charge_density_norm_squared
    )

    return charge_density_normalized

