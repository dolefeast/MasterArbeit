from math_objects.system import System
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
import scipy as sp
import numpy as np

def calculate_N_solutions(electric_potential: function, 
                          scalar_mass: float, 
                          scalar_charge: float, 
                          N: int,
                          omega_in: float,
                          omega_end: float,
                          scalar_value: float=1, 
                          tolerance: float=1e-2, 
                          n_points: int=300, 
                          scalar_name='phi'):
    """Calculates N  KG solutions associated to a certain external classical field
    Parameters
        electric_potential: function, Can be a scalar.
            If a function, electric_potential is the A0(z). If a float, 
            electric_potential is the value lambda in A0(z) = - lambda (z-1/2)
        scalar_mass: float, 
            The mass of the scalar field
        scalar_charge: float, 
            The charge of the scalar field
        N: int,
            Amount of desired solutions
        omega_in: float,
            initial omega for the sweep of eigenvalue
        omega_end: float,
            final omega for the sweep of eigenvalue
        scalar_value: float=1, 
            The initial guess of the scalar field. For solve_bvp purposes
        tolerance: float=1e-2, 
            The tolerance to check wether the solution converged
        n_points: int=300, 
            Number of nodes in the solution
        scalar_name='phi'
            Name of the scalar field. For representation purposes
    Returns:
        [scipy.integrate.solve_bvp.solution] of len() = N
    
    """
    system = System(field_strength=electric_potential,
            scalar_mass=scalar_mass,
            scalar_charge=scalar_charge,
            n_points=n_points,
            scalar_name=scalar_name,
            scalar_values=scalar_name)

    solution_array = []


    