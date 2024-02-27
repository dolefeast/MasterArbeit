import numpy as np
import scipy as sp
from math_objects.savitzky_golay import savitzky_golay


def calculate_total_charge_density(self, eigenstate_array, filter_window, filter_order):
    """
    Calculates the charge density in the system given a certain family of solutions
    Parameters:
        eigenstate_array: [scipy.integrate.solve_bvp solution objects]
    Returns:
        charge_density: np.array([float]) The total charge density at every node
    """
    e = self.phi.charge

    def total_charge_density(z):
        total_charge_density = 0
        cutoff_N = len(eigenstate_array) // 2
        for i, solution in enumerate(eigenstate_array):
            n = np.abs(i - cutoff_N)
            # we are adding and substracting the series -2eλ lim(t -> 0) sum z(1-z) sin(pi n z) cos (pi n z) exp(i pi n(t + i epsilon)) = -eλ z(1-z)/2 cot(pi z)
            # This way we can make the cutoff earlier, since each of the terms in the series will go to 0 faster.
            adding_term = -1
            renormalization_term = (
                e
                * self.lambda_value
                * z
                * (1 - z)
                * np.sin(np.pi * np.abs(n) * z)
                * np.cos(np.pi * n * z)
            )
            renormalization_term = savitzky_golay(renormalization_term, filter_window, filter_order)
            # print(f"Shape of smoothed out renormalization_term = {np.shape(renormalization_term)}")
            # print(f"Shape of calculate_charge_density = {np.shape(calculate_charge_density(solution)(z))}")
            total_charge_density = (
                total_charge_density
                + self.calculate_charge_density(solution)(z)
                + adding_term * renormalization_term
            )

        massless_charge_density = (
            e * self.lambda_value * z * (1 - z) / 2 / np.tan(np.pi * z)
        )
        potential_term = (
            1 / np.pi * e**2 * self.A0(z)
        )  # Extra term coming from point-split renormalization. Accounts for gauge invariance
        return (
            -total_charge_density
            + adding_term * massless_charge_density
            + potential_term
        )

    return total_charge_density
