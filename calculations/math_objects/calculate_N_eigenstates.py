from math_objects.system import System
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
import scipy as sp
import numpy as np


def calculate_N_eigenstates(
    electric_potential,
    scalar_mass: float,
    scalar_charge: float,
    N: int,
    omega_in: float,
    omega_end: float,
    boundary_conditions=None,
    scalar_value: float = 1,
    tolerance: float = 1e-2,
    n_points: int = 300,
    scalar_name: str = "phi",
    max_nodes: int = 3000,
):
    """Calculates N  KG eigenstates associated to a certain external classical field
    Parameters
        electric_potential: function, Can be a scalar.
            If a function, electric_potential is the A0(z).
            If scalar, electric_potential is the value lambda in A0(z) = - lambda (z-1/2)
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
        boundary_conditions: float=None,
            boundary conditions tos solve the differential equation.
            If not given (None), they are assumed to be dirichlet (y(0) = y(1) = 0)
        scalar_value: float=1,
            The initial guess of the scalar field. For solve_bvp purposes
        tolerance: float=1e-2,
            The tolerance to check wether the solution converged
        n_points: int=300,
            Number of nodes in the solution
        scalar_name: str='phi'
            Name of the scalar field. For representation purposes
        max_nodes: int=3000
            Maximum number of nodes for a given solution
    Returns:
        [scipy.integrate.solve_bvp.solution] of len() = N

    """
    if omega_in > 0:
        print(
            "Warning: The starting value for the omega guesses is positive, and the solution needs negative frequency solutions"
        )

    system = System(
        field_strength=electric_potential,
        scalar_mass=scalar_mass,
        scalar_charge=scalar_charge,
        n_points=n_points,
        scalar_name=scalar_name,
        scalar_values=scalar_name,
    )

    if boundary_conditions is None:
        boundary_conditions = system.dirichlet_boundary_conditions

    solution_array = []

    for i, omega_guess in enumerate(np.linspace(omega_in, omega_end, 500)):
        solution = sp.integrate.solve_bvp(
            system.differential_equation,
            boundary_conditions,
            system.z,
            (system.phi.value, system.phi.gradient.value),
            p=(omega_guess,),
            verbose=0,
            max_nodes=max_nodes,
            tol=TOL,
        )

        solution_array.append(solution)

    return solution_array
