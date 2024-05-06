from mpmath import quad, sqrt 
from scipy.interpolate import CubicSpline

def eigenstate_norm(self, eigenvalue, eigenstate):
    assert isinstance(eigenvalue, float), "eigenvalue must be a float"
    try:
        # Is eigenstate a solve_bvp solution?
        eigenstate.sol
        # It is, but I need eigenstate to be only the field, not its derivative.
        def eigenstate_callable(z):
            return eigenstate.sol(z)[0]

    except AttributeError:
        # It was not, it must be an array
        eigenstate_callable = CubicSpline(
                self.z, 
                eigenstate,
                )
        
    def rho_without_normalization(z):
        return abs(
            (eigenvalue - self.e * self.A0_field(z)) 
            * abs(eigenstate_callable(z))**2
            )

    norm_squared = quad(
        rho_without_normalization,
        [0, 1],
        )
    return sqrt(norm_squared)

def normalize_eigenstates(self):
    for i, (eigenvalue, eigenstate, eigenstate_gradient) in enumerate(
            zip(
                self.eigenvalue_array,
                self.eigenstate_array,
                self.eigenstate_gradient_array
                    )
                ):
        norm = self.eigenstate_norm(eigenvalue, eigenstate)
        # Convert it back to array format
        try:
            eigenstate = eigenstate.sol(self.z)[0]
        # It already was in array format
        except AttributeError:
            pass

        self.eigenstate_array[i] = eigenstate/norm
        self.eigenstate_gradient_array[i] = eigenstate_gradient/norm
