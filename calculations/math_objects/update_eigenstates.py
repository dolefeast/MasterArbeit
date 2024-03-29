import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from math_objects.get_eigenvalues import get_eigenvalues
from math_objects.fields import Vector_Potential
from math_objects.modify_A0 import modify_A0

from scripts.bao_filtering import bao_filtering
from scripts.float_to_str import float_to_str

def root_mean_square(x):
    return np.sqrt(
            np.mean(
                np.square(
                    x
                    )
                )
            )

def update_eigenstates(
        self,
        renormalization: bool=False,
        smoothing: bool=True,
        filter_method: callable=None,
        filter_parameters: tuple=None,
        save_results: bool=True,
        sig_digs: int=2,
        ):
    """
    Given the state of the system, returns the corresponding eigenstates off it
    Parameters:
        renormalization: bool=False. Whether we want 'renormalization' to be applied
        smoothing: bool=False. Whether we want the output of rho to be filtered
        filter_method: callable=None. call signature filter_method(x, signal, *filter_parameters) The filtering method to be used on the charge density. 
        filter_parameters: tuple=None. The paraeters to be passed to the filter_method
    """

    # Calculate the eigenstates of the actual state of the system,
    # given present A0 value
    eigenstate_array = self.calculate_eigenstates()

    # Calculate corresponding charge density
    total_rho_array = self.calculate_total_rho(
                renormalization=renormalization,
                )(self.z) # calculate_total_rho returns a function

    if smoothing:
        if filter_method is None:
            if filter_parameters is None:
                filter_parameters = (0.06,)
            print("No filtering method was given, but smoothing=True. Filtering using bao_smoothing...")
            z, total_rho_array = bao_filtering(
                    self.z, 
                    total_rho_array,
                    size=filter_parameters[0]
                    )[1]
        else:
            z, total_rho_array = filter_method(
                    self.z,
                    total_rho_array,
                    *filter_parameters
                    )
    else:
        z = np.linspace(0, 1, len(total_rho_array))

    # Since the filtering sometimes changes the mesh array
    #z, total_rho_array = total_rho_array
    
    # The mid-step calculations are allowed to change the shape of the arrays
    # Calculate the A0 corresponding to the new charge distribution
    A0_perturbation = modify_A0(z, total_rho_array)

    self.rho_array = total_rho_array

    self.A0_perturbation = A0_perturbation.sol(self.z)[0] - A0_perturbation.sol(0.5)[0]# Back to original shape and shifted to go to 0 at z=1/2

    self.A0_value = (
            - self.lambda_value * (self.z - 1/2) 
            + self.A0_perturbation
            )

    self.A0_field = Vector_Potential(
            value = (
        -self.lambda_value * (self.z - 1/2)
         + 1 * (
        + A0_perturbation.sol(self.z)[0] # This adds nothing since both arrays point to the same memory allocation. Each change in A0 field will change the value at self.A0_field.value
        - A0_perturbation.sol(0.5)[0] # To keep the problem antisymmetric wrt z=1/2
             )
            ),
            n_points=self.n_points
            )


def update_eigenstates_iteration(
    self,
    n_iterations:int=3,
    renormalization: bool=False,
    smoothing: bool=True,
    filter_method: callable=None,
    filter_parameters: tuple=None,
    save_results: bool=True,
    sig_digs: int=2,
    plot: bool=False,
    axis=None,
    ):
    """
    Updates the eigenstates iteratively
    Parameters:
    n_iterations;int =3, how many iterations of the routine
    renormalization: bool=False, whether to apply the trick
    smoothing: bool=True, if the total_rho_array is to be smoothed
    filter_method: callable=None, the filtering method to be used if smoothing=True
    filter_parameters: tuple=None, the filter parameters to pass to the filter_method
    save_results: bool=True, if the results are to be saved
    """

    for i in range(n_iterations):
        self.update_eigenstates(
        renormalization,
        smoothing,
        filter_method,
        filter_parameters,
        save_results=save_results,
        sig_digs=sig_digs,
                )
        if plot:
            if not axis is None:
                alpha = 0.3 + 0.7 * ((i+1)/n_iterations) ** 3
                axis.plot(
                        self.z,
                        self.rho_array,
                        'b',
                        alpha=alpha
                        )
                axis.plot(
                        self.z,
                        self.A0_perturbation,
                        'r',
                        alpha=alpha
                        )


def update_eigenstates_until_convergence(
    self,
    tol:float=1e-2,
    renormalization: bool=False,
    smoothing: bool=True,
    filter_method: callable=None,
    filter_parameters: tuple=None,
    save_results: bool=True,
    sig_digs: int=2,
    plot: bool=False,
    axis=None,
    ):
    """
    Updates the eigenstates iteratively until convergence is reached
    Parameters:
    tol:float=1e-2, tolerance for which convergence is reached (self.rho_array - previous_rho)/previous_rho
    renormalization: bool=False, whether to apply the trick
    smoothing: bool=True, if the total_rho_array is to be smoothed
    filter_method: callable=None, the filtering method to be used if smoothing=True
    filter_parameters: tuple=None, the filter parameters to pass to the filter_method
    save_results: bool=True, if the results are to be saved
    """

    # update the rho to create rho_array

    self.update_eigenstates(
        renormalization,
        smoothing,
        filter_method,
        filter_parameters,
        save_results=save_results,
        sig_digs=sig_digs,
            )

    previous_rho = np.copy(self.rho_array)

    self.update_eigenstates(
        renormalization,
        smoothing,
        filter_method,
        filter_parameters,
        save_results=save_results,
        sig_digs=sig_digs,
            )
    
    r = root_mean_square(
            (self.rho_array - previous_rho)
            /(1+previous_rho)
            )

    count = 1
    while r>tol:
        previous_rho = np.copy(self.rho_array)
        if plot:
            if not axis is None:
                axis.plot(
                        self.z,
                        self.rho_array,
                        'b',
                #        alpha=alpha
                        )
                axis.plot(
                        self.z,
                        self.A0_perturbation,
                        'r',
                #        alpha=alpha
                        )
        self.update_eigenstates(
            renormalization,
            smoothing,
            filter_method,
            filter_parameters,
            save_results=save_results,
            sig_digs=sig_digs,
        )
        count += 1
        r = root_mean_square(
                    (self.rho_array - previous_rho)
                    /(1+previous_rho)
                )

    if save_results and not self.broken:
        self.save_solutions()

    if plot:
        if not axis is None:
            axis.plot(
                    self.z,
                    self.rho_array,
                    'b',
            #        alpha=alpha
                    )
            axis.plot(
                    self.z,
                    self.A0_perturbation,
                    'r',
            #        alpha=alpha
                    )

    print(f"Convergence was reached after {count} iterations")

