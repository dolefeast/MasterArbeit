import numpy as np
from math_objects import Vacuum_Solution

from scripts.read_files import read_files
from scripts.extrapolation import linear_extrapolation
from scripts.get_lambda_index import get_lambda_index

def init_axes(
        plot,
        ):
    # The idea is: At different times I will want 
    # to plot different things. I want to create
    # different fig environments to plot each of these
    # things
    if plot:
        global plt
        import matplotlib.pyplot as plt

        global fig, ax_fields, ax_omegas
        fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))

        ax_fields.set_xlabel("z")
        ax_fields.set_ylabel(r"$\rho(z)$")

        ax_omegas.set_xlabel("$\lambda$")
        ax_omegas.set_ylabel(r"$max(A_0)$")
    else:
        ax_fields = None
        ax_omegas = None

def main(
        m,
        lambda_step=None,
        lambda_iterations=0,
        lambda_min:float=None,
        plot_max_A0:bool=False,
        N_mode_cutoff:int=49,
        n_points:int=400,
        n_iterations=None,
        e:float=1,
        bcs:str='dirichlet',
        burnout:int=0,
        plot=True,
        save_results=True,
        ):


    init_axes(plot=plot)

    data = read_files(m=m)
    (
        eigenvalue_array,
        eigenstate_array,
        eigenstate_gradient_array,
        A0_perturbation_file_array,
        rho_file_array,
        lambda_value_array
    ) = data.values()
    if lambda_min is not None:  
        # To validate that lambda_min was a float
        try: 
            float(lambda_min)
            # To get the lambda value that is just smaller than the given lambda_min
            burnout = len(lambda_value_array) - get_lambda_index(lambda_value_array, lambda_min)
        except ValueError:
            raise ValueError('Argument min_lambda_was given, but it was not a scalar')

    # Burnout because studying for lambda too close to the 
    # maximum lambda may not give enough results to draw 
    # conclusions
    eigenvalue_array=eigenvalue_array[:-1-burnout]
    A0_perturbation_file_array=A0_perturbation_file_array[:-1-burnout]
    rho_file_array=rho_file_array[:-1-burnout]
    lambda_value_array=lambda_value_array[:-1-burnout]
    
    max_A0_array = [
            max(A0_perturbation_file)
            for A0_perturbation_file in A0_perturbation_file_array
            ]
    
    print()
    system = Vacuum_Solution(
        lambda_value=max(lambda_value_array),
        scalar_name="phi",
        n_points=n_points,
        N_mode_cutoff=N_mode_cutoff + 1,
        e=e,
        m=m,
        float_tol=1e-1,
        read_solutions=True,
        bcs=bcs,
    )

    next_lambda_value, next_A0_factor = linear_extrapolation(
            lambda_value_array,
            max_A0_array,
            x_step=lambda_step
            )
    print("\n----- Starting calculations -----")

    lambda_count = 0
    try:
        while True:  
            print(f"\nIncreasing lambda to {next_lambda_value}")
            # Update the field
            previous_lambda_value = system.lambda_value
            ax_omegas.plot(
                    system.lambda_value,
                    max(
                        system.A0_induced
                       # system.A0_field.value 
                       # + system.lambda_value * (system.z - 1/2)),
                       ),
                    '^k'
                    )

            next_lambda_value, next_A0_factor = linear_extrapolation(
                    lambda_value_array,
                    max_A0_array,
                    x_step=lambda_step
                    )
            #system.A0_induced = system.A0_induced * next_lambda_value/system.lambda_value

            system.A0_induced *= 1.5*next_A0_factor/max(system.A0_induced)
            ax_omegas.plot(
                    system.lambda_value,
                    max(
                        system.A0_induced
#                        system.A0_field.value 
#                        + system.lambda_value 
#                        * (system.z - 1/2)
                        ),
                    'xr'
                    )

            system.lambda_value = next_lambda_value
            # Calculate the solutions for the new external E
            system.update_eigenstates_script(
                    n_iterations=n_iterations,
                    plot=plot,
                    axis=ax_fields,
                    save_results=save_results,
                )
            if system.broken:
                break
            if plot:
                ax_omegas.plot(
                        previous_lambda_value,
                        max(
                            system.A0_field.value 
                            + system.lambda_value 
                            * (system.z - 1/2)
                            ),
                        'go'
                        )

            max_A0_array.append(max(system.A0_induced))
            lambda_value_array.append(next_lambda_value)
            next_lambda_value, next_A0_factor = linear_extrapolation(
                    lambda_value_array,
                    max_A0_array/max_A0_array[0],
                    x_step=lambda_step
                    )
            if lambda_iterations is not None:
                if lambda_count>lambda_iterations:
                    break
                lambda_count+=1
    except KeyboardInterrupt:
        pass

    if plot:
        ax_omegas.plot([], [], '^k', label='Unmodified previous max $A_0$')
        ax_omegas.plot([], [], 'xr', label='Modified previous solution')
        ax_omegas.plot([], [], 'go', label='Actual solution')

        fig.suptitle(f"$m={m}$")
        ax_fields.set_xlabel("z")
        ax_fields.set_ylabel(r"$\rho(z)$")

        ax_omegas.set_xlabel("$\lambda$")
        ax_omegas.set_ylabel(r"$max(A_0)$")
        ax_omegas.legend(loc=r"best")

        fig.tight_layout()
        plt.show()


if __name__ == "__main__":
    N_mode_cutoff = 49 # Resulting rho(z) has freq of 1/(N+1)
    n_points = (N_mode_cutoff+1)*8
    n_iterations = 3

    lambda_iterations = None
    main(
            m=3.,
        lambda_min=20.5,
        lambda_step=0.7,
        n_iterations=n_iterations,
        lambda_iterations=lambda_iterations,
        N_mode_cutoff=N_mode_cutoff, 
        n_points=n_points,
        bcs='dirichlet',
        burnout=20,
        plot=True,
        save_results=False,
            )
