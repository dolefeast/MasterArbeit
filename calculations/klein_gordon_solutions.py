from math_objects import Vacuum_Solution
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize

from scripts.plotting import plot_each_eigenstate, plot_different_window_filter, scatter_omegas, plot_from_0_to_1
import scripts.filter_scripts as filter_scripts
import scripts.filters as filters
from scripts.bao_filtering import bao_filtering

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
from pathlib import Path

def init_main(
        N_mode_cutoff,
        lambda_value,
        TOL,
        N_POINTS,
        m,
        e,
        smoothing=False,
        plot=True,
        read_solutions=True,
        ):
    # Initialize physical system
    system = Vacuum_Solution(
        lambda_value=lambda_value,
        scalar_name="phi",
        n_points=N_POINTS,
        N_mode_cutoff=N_mode_cutoff+1,
        e=e,
        m=m,
        float_tol=1e-1,
        read_solutions=read_solutions,
    )

    if plot: # Just to fold this
        global plt
        import matplotlib.pyplot as plt
        global fig, ax_fields, ax_omegas
        fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))
        fig.suptitle(f"$m={m}$")

        ax_fields.set_xlabel('z')
        ax_fields.set_ylabel(r'$\rho(z)$')

        ax_omegas.set_xlabel('n')
        ax_omegas.set_ylabel(r'$\omega_n$')

        ax_omegas.plot([],[], 'go', label='positive eigenvalue')
        ax_omegas.plot([],[], 'rx', label='negative eigenvalue')
#        ax_fields.legend(loc='best')
#        ax_omegas.legend(loc='best')
        plt.tight_layout()
    else: 
        ax_fields = None
        ax_omegas = None

        
    return system

def main(
        N_mode_cutoff,
        lambda_min,
        lambda_max,
        lambda_div,
        TOL,
        N_POINTS,
        m,
        e,
        n_iterations,
        save_results=True,
        smoothing=True,
        plot=True,
        read_solutions=True,
    ):

    system = init_main(
        N_mode_cutoff,
        lambda_min,
        TOL,
        N_POINTS,
        m,
        e,
        plot=plot,
        read_solutions=read_solutions
        )

    filter_method = filter_scripts.extend_and_filter
    filter_parameters = (
            filters.convolve_twice,
            tuple((1, )),
            9/N_POINTS,
            )

    for i, iterating_lambda in enumerate(
            np.linspace(
                lambda_min,
                lambda_max,
                lambda_div
            )
            ):
        print(20*'=')
        print(f'iterating_lambda = {iterating_lambda}')
        system = Vacuum_Solution(
                lambda_value=iterating_lambda,
                A0_perturbation= (
                    (
                        system.A0_field.value 
                        + system.lambda_value * (system.z - 1/2)
                        )
                    *iterating_lambda
                    /system.lambda_value
                    ),
                m=m,
                e=e,
                n_points=N_POINTS,
                eigenvalue_array=system.eigenvalue_array,
                eigenstate_array=system.eigenstate_array,
                eigenstate_gradient_array=system.eigenstate_gradient_array,
                float_tol=1e-2,
                read_solutions=True,
                )

                
        system.update_eigenstates_until_convergence(
                #n_iterations=2,
                tol=1e-1,
                smoothing=smoothing,
                filter_method=filter_method,
                filter_parameters=filter_parameters,
                save_results=save_results,
                plot=plot,
                axis=ax_fields,
                #  filter_method=filtering.double_filtering,
                #  filter_parameters=(150,)
                )
        if system.broken:
            print("The calculation broke in the previous iteration. Breaking iteration...")
            break
        x_density, y_density = plot_from_0_to_1(system.rho_array)

        x_field, y_field = plot_from_0_to_1(
                system.A0_field(system.z) 
                + iterating_lambda * (system.z-1/2) # minus the base value
            ) # to see perturbation

        # The higher the exponent, the more the first iterations will be'muted'
        alpha = 0.3 + 0.7*((i+1)/lambda_div)**5


    if plot:
        scatter_omegas(system.eigenvalue_array, ax_omegas, m)
#        ax_fields.legend(loc='best')
        plt.show()
    

if __name__ == "__main__":
    from scripts.input_list import input_list

    N_mode_cutoff = 49 # Resulting rho(z) has freq of 1/(N+1)
    N_POINTS = (N_mode_cutoff+1)*8

    TOL = 1e-2
    e = 1
    n_iterations = None

    lambda_min = 20.861
    lambda_max = 25
    lambda_div = 25
    for m in np.linspace(3, 10, 1):
        system = main(
            N_mode_cutoff,
            lambda_min,
            lambda_max,
            lambda_div,
            TOL,
            N_POINTS,
            m,
            e,
            n_iterations,
            smoothing=True,
            save_results=True,
            read_solutions=True,
            plot=True,
        )

