from math_objects import Vacuum_Solution
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.savitzky_golay import savitzky_golay

from scripts.plotting import plot_each_eigenstate, plot_different_window_filter, scatter_omegas, plot_from_0_to_1
import scripts.filtering as filtering
from scripts.bao_filtering import bao_filtering

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def init_main(
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

    if True: # Just to fold this
        global fig, ax_fields, ax_omegas
        fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))
        fig.suptitle(f"$\lambda={lambda_value}, m={m}$")

        ax_fields.set_xlabel('z')
        ax_fields.set_ylabel(r'$\rho(z)$')

        ax_omegas.set_xlabel('n')
        ax_omegas.set_ylabel(r'$\omega_n$')

        ax_omegas.plot([],[], 'go', label='positive eigenvalue')
        ax_omegas.plot([],[], 'rx', label='negative eigenvalue')
        ax_fields.legend(loc='best')
        ax_omegas.legend(loc='best')
        plt.tight_layout()
        
    return system

def main(
        N_mode_cutoff,
        lambda_value,
        TOL,
        N_POINTS,
        m,
        e,
        smoothing=False,
    ):

    system = init_main(
        N_mode_cutoff,
        lambda_value,
        TOL,
        N_POINTS,
        m,
        e,
        )

    n_interations = 1 # N of iterations to update the electric potential

    for index, i in enumerate(np.linspace(0.3, 1, n_interations)):
        print(f'Iteration number {index+1}')
        system.update_eigenstates(
                smoothing=True,
                filter_parameters=(0.12,)
                #  filtering_method=filtering.double_filtering,
                #  filter_parameters=(150,)
                )
        x_density, y_density = plot_from_0_to_1(system.charge_density_array)
        x_field, y_field = plot_from_0_to_1(system.A0_perturbation)
        ax_fields.plot(x_density, y_density, 'b', alpha=0.3 + 0.3*i)
        # ax_fields.plot(x_field, y_field, 'r', alpha=0.3 + 0.3*i)

    # to_csv = np.asarray((x_density, y_density))
    # np.savetxt(f'./saved_solutions/lambda_{lambda_value}_mass_{m}.csv', to_csv, delimiter=",")
    
    scatter_omegas(system.eigenvalue_array, ax_omegas, m)
    plt.show()

if __name__ == "__main__":

    N_mode_cutoff = 50
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
