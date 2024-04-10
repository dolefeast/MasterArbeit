from math_objects.unique_floats import float_in_array
from scipy import integrate
import numpy as np

import mpmath

def normalize_eigenstates(self, true_eigenstate):
    """Given a PPoly sol object, normalize it"""
    true_eigenvalue = true_eigenstate.p[0]
    rho_without_normalization = lambda z: (
        (true_eigenvalue - self.e * self.A0_field(z))
        * np.abs(true_eigenstate.sol(z)[0]) ** 2
    )

    norm_squared = mpmath.quad(rho_without_normalization, [0, 1])

    norm = np.sqrt(np.abs(norm_squared))

    eigenstate_normalized = true_eigenstate.sol(self.z)[0] / norm
    eigenstate_gradient_normalized = true_eigenstate.sol(self.z)[1] / norm
    rho_n_normalized = (
        np.sign(true_eigenvalue) * rho_without_normalization(self.z) / norm_squared
    )
    return eigenstate_normalized, eigenstate_gradient_normalized, rho_n_normalized


def calculate_eigenstates(
    self,
    float_tol: float = 1e-2,
    bvp_tol: float = 5e-5,
    max_nodes: int = 3000,
    verbose: int = 0,
):
    """Calculates KG eigenstate_array of len(eigenstate_guess) associated to a certain external classical field.
    Parameters:
        float_tol: float=1e-2. The tolerance for which two floats are considered the same
        bvp_tol: float=1e-2. The tolerance for which the KG are submitted to
        max_nodes: int=3000. The maximal amount of nodes solve_bvp is allowed to have
        verbose: int=0. Can take values in [0,1,2]. verbosity in increasing order
    """
    solution_array = []
    solution_gradient_array = []
    true_eigenvalue_array = []
    rho_n_array = []

    repeated_eigenvalue_count = 0

    for eigenstate_guess, eigenstate_gradient_guess, eigenvalue_guess in zip(
        self.eigenstate_array,
        self.eigenstate_gradient_array,
        self.eigenvalue_array
    ):
        # This never converges, it is a non physical solution
        if eigenvalue_guess == 0.0:
            continue

        true_eigenstate = integrate.solve_bvp(
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
        # Check if eigenvalue is double counted, or if solve_bvp converged
        if float_in_array(true_eigenvalue, true_eigenvalue_array, tol=self.float_tol):
            # Discarding either repeated eigenvalues or non convergent solutions
            print(
                f"Found repeated eigenvalue: {true_eigenstate.p[0]} with\neigenvalue_guess={eigenvalue_guess}. Escaping iteration... \n"
            )
            repeated_eigenvalue_count += 1
            self.broken = 1
            break
        if not true_eigenstate.success and eigenvalue_guess != 0.0:
            if abs(eigenvalue_guess) < 0.4:
                print(
                    f"Warning: eigenvalue_guess={eigenvalue_guess} did not converge.\n\tRemoving the following eigenvalues (and their respective eigenstates) and starting again"
                        )

                n = len(self.eigenstate_array) 
                print("\t", self.eigenvalue_array.pop(n//2-1))
                print("\t", self.eigenvalue_array.pop(n//2-1))
                self.eigenstate_array.pop(n//2-1)
                self.eigenstate_array.pop(n//2-1)
                self.eigenstate_gradient_array.pop(n//2-1)
                self.eigenstate_gradient_array.pop(n//2-1)
                self.calculate_eigenstates()
            else:
                print(
                    f"Warning: Eigenvalue={eigenvalue_guess} did not converge.\n\tEscaping iteration... "
                )
                print("Status of the bvp solution:", true_eigenstate.message)
                self.broken = 1
                break
        elif len(solution_array) >= len(self.eigenstate_array)//2:
            pass


        (
                eigenstate_normalized,
                eigenstate_gradient_normalized,
                rho_n_normalized
                ) = normalize_eigenstates(self, true_eigenstate)

        solution_array.append(eigenstate_normalized)
        solution_gradient_array.append(eigenstate_gradient_normalized)
        true_eigenvalue_array.append(true_eigenvalue)
        rho_n_array.append(rho_n_normalized)

        if verbose:
            print(true_eigenstate)

        # true_eigenvalues_array.append(true_eigenvalue)

    self.eigenstate_array = solution_array
    self.eigenstate_gradient_array = solution_gradient_array
    self.eigenvalue_array = true_eigenvalue_array
    self.rho_n_array = rho_n_array

    # return solution_array
