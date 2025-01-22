from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from mpmath import quad
import numpy as np


def calculate_eigenstates(self):
    """
    Calculates the solutions to the klein_gordon equation by parametrizing the solution family by a parameter, and then looking for the parameter that verifies the boundary conditions
    """

    if self.bcs == "dirichlet":
        initial_values = (0, 1)
        bcs_index = 0 
    elif self.bcs == "neumann":
        initial_values = (1, 0)
        bcs_index = 1 

    eigenvalue_lower_bound = self.bisection_method_lower_bound(self.eigenvalues)
    eigenvalue_upper_bound = self.bisection_method_upper_bound(self.eigenvalues)

    parametrized_ODE = lambda omega: solve_ivp(
            lambda z, y: self.Klein_Gordon_equation(z, y, omega), (0,1), initial_values, dense_output=True
            ) # The first 0 1 is the range, the second one are the initial values

    # import matplotlib.pyplot as plt


    # omega_array = np.arange(0, 20, 0.1)
    # plt.plot(omega_array, [parametrized_ODE(omega).sol(1)[1] for omega in omega_array], )
    # plt.show()

    # exit()
    eigenvalues = [
            self.find_root_bisection(
                lambda omega: parametrized_ODE(omega).sol(1)[bcs_index], *omega_upper_lower
                ) 
            for omega_upper_lower in zip(
                eigenvalue_lower_bound,
                eigenvalue_upper_bound
                )
            ]

    # the eigenvalues should be antisymmetric i.e. omega_n = -omega_{-n}
    self.eigenvalues = [ (i-j) / 2 for i, j in zip(eigenvalues, eigenvalues[::-1]) ]
    self.eigenstates = [ parametrized_ODE(omega).sol(self.z)[bcs_index] for omega in eigenvalues ]

def normalize_eigenstates(self):
    """
    Normalizes the eigenstates 
    """
    # Warning, math

    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
                self.eigenvalues,
                self.eigenstates,
            )
            ):

        eigenstate = CubicSpline(self.z, eigenstate)

        def rho_n_without_normalizing(z):
            # Normalizing wrt the symplectic norm
            # the solutions need not be real.
            return (eigenvalue - self.A0(z)) * abs(eigenstate(z))**2

        # Calculate the norm
        norm_squared = abs(float(quad(rho_n_without_normalizing, [0, 1])))

        norm = self.sign(eigenvalue)*np.sqrt(norm_squared)

        self.eigenstates[n] /= norm
