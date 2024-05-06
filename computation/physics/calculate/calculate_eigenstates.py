import scipy.integrate as integrate

from utils.float_in_array import float_in_array

def calculate_eigenstates(
    self,
    float_tol:float=1e-2,
    bvp_tol:float=1e-5,
    max_nodes:int=2e5,
    verbose:int=0,
    ):
    """
    Calculates the eigenstate family 
    Parameters:
        bvp_tol:float=1e-2, the tolerance to which solutions are calculated
        float_tol:float=1e-5, the tolerance to which to floats are said to be different
        max_nodes:int=2e5, max_nodes for solve_bvp
        verbose:int=False,	controls verbosity of solve_bvp
    """

    # true_ distinction means they are solutions
    # whereas self.*_array are just fiducial values (in this case)
    true_eigenstate_array = []
    true_eigenstate_gradient_array = []
    true_eigenvalue_array = []

    for i, (
            eigenstate_guess,
            eigenstate_gradient_guess,
            eigenvalue_guess
            ) in enumerate(
            zip(
                self.eigenstate_array,
                self.eigenstate_gradient_array,
                self.eigenvalue_array
                )
                    ):

        # Non physical
        if eigenvalue_guess == 0:
            continue

        true_eigenstate = integrate.solve_bvp(
            self.Klein_Gordon,
            self.boundary_conditions,
            self.z,
            (eigenstate_guess, eigenstate_gradient_guess),
            p=(eigenvalue_guess,),
            verbose=0,
            max_nodes=max_nodes,
            tol=bvp_tol,
        )

        true_eigenvalue = true_eigenstate.p[0]

        # First check if this eigenvalue has already been solved for
        if float_in_array(true_eigenvalue,
                true_eigenvalue_array,
                tol=self.float_tol
                ): 
            print(
                f"Found repeated eigenvalue: {true_eigenstate.p[0]} with\neigenvalue_guess={eigenvalue_guess}. Escaping iteration... \n"
                )
            self.broken = 1
            break
        
        # If the calculation broke down stop
        if not true_eigenstate.success:
            print(
            f"Warning: Eigenvalue={eigenvalue_guess} did not converge.\n\tEscaping iteration... "
            )

            self.broken = 1
            break

        # If it did not, save the solution
        true_eigenvalue_array.append(true_eigenvalue)
        # Be careful here. In the true eigenstate_array it saves the solve_bvp
        # object, but in the true_eigenstate_gradient_array it saves the array.
        # The true_eigenstate_array MUST be later modified.
        true_eigenstate_array.append(true_eigenstate) 
        true_eigenstate_gradient_array.append(true_eigenstate.sol(self.z)[1])

    self.eigenvalue_array = true_eigenvalue_array
    self.eigenstate_array = true_eigenstate_array
    self.eigenstate_gradient_array = true_eigenstate_gradient_array
