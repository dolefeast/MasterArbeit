import numpy as np
import scipy as sp

def calculate_N_eigenstates(
    self,
    omega_in: float,
    omega_end: float,
    N: int,
    boundary_conditions=None,
    tolerance: float = 1e-2,
    n_points: int = 300,
    max_nodes: int = 3000,
):
    """Calculates N (not necessarily normalized) KG eigenstates associated to a certain external classical field
    Parameters
        N: int,
            Amount of desired eigenstates.
        omega_in: float,
            initial omega for the sweep of eigenvalue
        omega_end: float,
            final omega for the sweep of eigenvalue
        boundary_conditions: float=None,
            boundary conditions tos solve the differential equation.
            If not given (None), they are assumed to be dirichlet (y(0) = y(1) = 0)
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

    if boundary_conditions is None:
        boundary_conditions = self.dirichlet_boundary_conditions

    solution_array = []
    #print(f'omega_in={omega_in}, omega_end={omega_end}, N={N}')

    for i, omega_guess in enumerate(np.linspace(omega_in, omega_end, N)):
        eigenfunction_guess = np.sin(
            omega_guess * self.z
        )  # Guessing the solution of the ODE
        guess_derivative = np.diff(eigenfunction_guess)
        guess_derivative = np.append(
            guess_derivative, guess_derivative[-1]
        )  # since calculating derivatives substracts the last point

        solution = sp.integrate.solve_bvp(
            self.differential_equation,
            boundary_conditions,
            self.z,
            (eigenfunction_guess, guess_derivative),
            p=(omega_guess,),
            verbose=0,
            max_nodes=max_nodes,
            tol=tolerance,
        )

        if solution.success:  # If it converged, add it to the solutions array
            solution_array.append(solution)

    return solution_array
