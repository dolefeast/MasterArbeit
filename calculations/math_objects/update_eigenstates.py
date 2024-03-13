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
        filter_parameters: tuple=None
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
            eigenstate_array, 
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
    self.eigenstate_array = [eigenstate.sol(self.z)[0] for eigenstate in eigenstate_array]
    self.eigenstate_gradient_array = [eigenstate.sol(self.z)[1] for eigenstate in eigenstate_array]
    self.eigenvalue_array = get_eigenvalues(eigenstate_array)

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

