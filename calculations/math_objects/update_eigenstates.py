import numpy as np
import scipy as sp

from math_objects.get_eigenvalues import get_eigenvalues
from math_objects.fields import Vector_Potential
from math_objects.modify_A0 import modify_A0

from scripts.bao_filtering import bao_filtering

def update_eigenstates(
        self,
        renormalization: bool=False,
        smoothing: bool=True,
        filter_method: callable=None,
        filter_parameters: tuple=None,
        save_results: bool=True,
        ):
    """
    Given the state of the system, returns the corresponding eigenstates off it
    Parameters:
        renormalization: bool=False. Whether we want 'renormalization' to be applied
        smoothing: bool=False. Whether we want the output of charge_density to be filtered
        filter_method: callable=None. call signature filter_method(x, signal, *filter_parameters) The filtering method to be used on the charge density. 
        filter_parameters: tuple=None. The paraeters to be passed to the filter_method
    """

    # Calculate the eigenstates of the actual state of the system,
    # given present A0 value
    eigenstate_array = self.calculate_eigenstates()

    # Calculate corresponding charge density
    total_charge_density_array = self.calculate_total_charge_density(
                renormalization=renormalization,
                )(self.z) # calculate_total_charge_density returns a function

    if smoothing:
        if filter_method is None:
            if filter_parameters is None:
                filter_parameters = (0.06,)
            print("No filtering method was given, but smoothing=True. Filtering using bao_smoothing...")
            total_charge_density_array = bao_filtering(
                    self.z, 
                    total_charge_density_array,
                    size=filter_parameters[0]
                    )[1]
        else:
            print("Filtering with the given filter_method")
            total_charge_density_array = filter_method(
                    self.z,
                    total_charge_density_array,
                    *filter_parameters
                    )
    # Since the filtering sometimes changes the mesh array
    z, total_charge_density_array = total_charge_density_array
    
    # The mid-step calculations are allowed to change the shape of the arrays
    # Calculate the A0 corresponding to the new charge distribution
    A0_perturbation = modify_A0(z, total_charge_density_array)

    self.charge_density_array = total_charge_density_array

    self.A0_perturbation = A0_perturbation.sol(self.z)[0] - A0_perturbation.sol(0.5)[0]# Back to original shape and upscaled it to go to 0

    self.A0_value = self.A0_base_value + self.A0_perturbation 


    self.A0_field = Vector_Potential(
            value = (
        -self.lambda_value * (self.z - 1/2)
        + A0_perturbation.sol(self.z)[0] # This adds nothing since both arrays point to the same memory allocation. Each change in A0 field will change the value at self.A0_field.value
        - A0_perturbation.sol(0.5)[0] # To keep the problem antisymmetric wrt z=1/2
            ),
            n_points=self.n_points
            )

    if save_results:
        lambda_string = str(float(self.lambda_value)).replace(".", "_")
        m_string = str(float(self.m)).replace(".", "_")
        file_id = f'lambda_{lambda_string}_mass_{m_string}.txt'
        to_csv = np.concatenate(((self.eigenvalue_array,), 
                            np.array(self.eigenstate_array).T)).T
        np.savetxt(f'saved_solutions/normalized_eigenstate/{file_id}', to_csv, delimiter= ",")

        to_csv = np.concatenate(((self.eigenvalue_array,), 
                            np.array(self.eigenstate_gradient_array).T)).T
        np.savetxt(f'saved_solutions/normalized_eigenstate_gradient/{file_id}', to_csv, delimiter= ",")


        to_csv = self.A0_field.value
        np.savetxt(f'saved_solutions/A0_field/{file_id}', to_csv, delimiter= ",")

        np.savetxt(f'saved_solutions/charge_density/{file_id}', to_csv, delimiter= ",")

def update_eigenstates_iteration(
    self,
    n_iterations:int=3,
    renormalization: bool=False,
    smoothing: bool=True,
    filter_method: callable=None,
    filter_parameters: tuple=None,
    save_results: bool=True,
    ):
    """
    Updates the eigenstates iteratively
    Parameters:
    n_iterations;int =3, how many iterations of the routine
    renormalization: bool=False, whether to apply the trick
    smoothing: bool=True, if the total_charge_density_array is to be smoothed
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
        save_results,
                )
