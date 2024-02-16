from math_objects import Vacuum_Solution
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.savitzky_golay import savitzky_golay

from scripts.plotting import plot_each_eigenstate, plot_different_window, scatter_omegas

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def moving_average(x, w): return np.convolve(x, np.ones(w), 'full') / w

def main(
        N_mode_cutoff,
        lambda_value,
        TOL,
        N_POINTS,
        m,
        e,
        smoothing=False,
    ):

    # Initialize physical system
    system = Vacuum_Solution(
        lambda_value=lambda_value,
        scalar_name="phi",
        n_points=N_POINTS,
        N_mode_cutoff=N_mode_cutoff+1,
        e=e,
        m=m,
        float_tol=1e-1
    )

    fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))
    fig.suptitle(f"$\lambda={lambda_value}, m={m}$")

    # system.phi.value = np.sin(system.z)

    eigenstates = []
    omega_solution_array = []

    eigenstate_array = system.calculate_eigenstates()

    z = system.z # Because I consistently make the mistake of only asking for z

    scatter_omegas(eigenstate_array, ax_omegas, m)
    total_charge_density = system.calculate_total_charge_density(
            eigenstate_array, 
                renormalization=False,
                smoothing=smoothing,
                )(z) # calculate_... returns a function

    ax_fields.plot(
            total_charge_density,
            label='Unfiltered charge density',
            alpha=0.4
            )

    plot_different_window(total_charge_density, ax_fields, moving_average)

    ax_fields.set_xlabel('z')
    ax_fields.set_ylabel(r'$\rho(z)$')

    ax_omegas.set_xlabel('n')
    ax_omegas.set_ylabel(r'$\omega_n$')

    ax_omegas.plot([],[], 'go', label='positive eigenvalue')
    ax_omegas.plot([],[], 'rx', label='negative eigenvalue')
    ax_fields.legend(loc='best')
    ax_omegas.legend(loc='best')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    N_mode_cutoff = 150
    lambda_value = 1
    TOL = 1e-2
    N_POINTS = 1000
    m = 5
    e = 1

    main(
        N_mode_cutoff,
        lambda_value,
        TOL,
        N_POINTS,
        m,
        e,
        smoothing=True,
    )
