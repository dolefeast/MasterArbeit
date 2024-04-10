from math_objects import Vacuum_Solution

from scripts.read_files import read_files, max_A0
from scripts.extrapolation import linear_extrapolation

def init_main(
        plot_max_A0,
        ):
    # The idea is: At different times I will want 
    # to plot different things. I want to create
    # different fig environments to plot each of these
    # things
    if plot_max_A0:
        import matplotlib.pyplot as plt
        fig, ax_max_A0 = plt.subplots()

def main(
        m,
        lambda_step=None,
        plot_max_A0:bool=False,
        N_mode_cutoff:int=49,
        n_points:int=400,
        n_iterations=None,
        e:float=1,
        bcs:str='dirichlet',
        burnout:int=0,
        plot=True,
        ):

    data = read_files(m=m)
    (
        eigenvalue_array,
        eigenstate_array,
        eigenstate_gradient_array,
        A0_perturbation_file_array,
        rho_file_array,
        lambda_value_array
    ) = data.values()


    eigenvalue_array=eigenvalue_array[:-1-burnout]
    A0_perturbation_file_array=A0_perturbation_file_array[:-1-burnout]
    rho_file_array=rho_file_array[:-1-burnout]
    lambda_value_array=lambda_value_array[:-1-burnout]
    
    max_A0_array = [
            max_A0(A0_perturbation_file)
            for A0_perturbation_file in A0_perturbation_file_array
            ]
    
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

#    system.update_eigenstates_script(
#            n_iterations=n_iterations,
#            plot=plot,
#            )


    next_lambda_value, next_A0_factor = linear_extrapolation(
            lambda_value_array,
            max_A0_array/max_A0_array[0],
            x_step=lambda_step
            )

    print(f"Increasing lambda to {next_lambda_value}")
    
    system = Vacuum_Solution(
        lambda_value=next_lambda_value,
        A0_induced=system.A0_induced * next_lambda_value/system.lambda_value,
        scalar_name="phi",
        n_points=n_points,
        N_mode_cutoff=N_mode_cutoff + 1,
        eigenvalue_array=system.eigenvalue_array,
        eigenstate_array=system.eigenstate_array,
        eigenstate_gradient_array=system.eigenstate_gradient_array,
        e=e,
        m=m,
        float_tol=1e-1,
        read_solutions=False,
        bcs=bcs,
    )


    system.update_eigenstates_script(
            n_iterations=n_iterations,
            plot=plot,
            )

if __name__ == "__main__":
    N_mode_cutoff = 49 # Resulting rho(z) has freq of 1/(N+1)
    n_points = (N_mode_cutoff+1)*8
    main(m=5.,
        lambda_step=None,
        N_mode_cutoff=N_mode_cutoff, 
            n_points=n_points,
            bcs='dirichlet',
            burnout=10,
            plot=True,
            )
