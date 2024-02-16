from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
import scipy as sp
import numpy as np

def calculate_eigenstates(self,
        float_tol: float=1e-2,
        bvp_tol: float=1e-2,
        max_nodes: int=3000,
        verbose: int=0,
        ):
    """Calculates len(self.eigenvalue_guess_array) KG eigenstate_array associated to a certain external classical field.
    Parameters:
        float_tol: float=1e-2. The tolerance for which two floats are considered the same
        bvp_tol: float=1e-2. The tolerance for which the KG are submitted to
        max_nodes: int=3000. The maximal amount of nodes solve_bvp is allowed to have
        verbose: int=0. Can take values in [0,1,2]. verbosity in increasing order
        """
    solution_array = []
    true_eigenvalues_array = []

    repeated_eigenvalue_count = 0

    for eigenstate_guess, eigenstate_gradient_guess, eigenvalue_guess in zip(
            self.eigenstate_array,
            self.eigenstate_gradient_array,
            self.eigenvalue_guess_array
            ):

        true_eigenstate = sp.integrate.solve_bvp(
            self.Klein_Gordon,
            self.boundary_conditions,
            self.z,
            (eigenstate_guess, eigenstate_gradient_guess),
            p=(eigenvalue_guess,),
            verbose=verbose,
            max_nodes=max_nodes,
            tol=bvp_tol,
        )

        true_eigenvalue = true_eigenstate.p[0]
        # Check if eigenvalue is double counted, or if it converged.
        if float_in_array(true_eigenvalue, true_eigenvalues_array, tol=self.float_tol):
            # Discarding either repeated eigenvalues or non convergent solutions
            print(f'Found repeated eigenvalue: {true_eigenstate.p[0]} with\neigenvalue_guess={eigenvalue_guess}\n')
            repeated_eigenvalue_count += 1
            continue
        if not true_eigenstate.success:
            print(f'Eigenvalue={eigenvalue_guess} did not converge')
            continue

        if verbose:
            print(true_eigenstate)

        solution_array.append(true_eigenstate)
        true_eigenvalues_array.append(true_eigenvalue)

    print(f"Found and removed {repeated_eigenvalue_count} repeated eigenvalues!")

    return solution_array


