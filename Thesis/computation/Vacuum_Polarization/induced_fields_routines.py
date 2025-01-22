import numpy as np
from scipy.interpolate import CubicSpline


def calculate_rho(self, ):
    """
    Calculates the vacuum polarization AFTER normalizing the eigenstates, but before filtering
    """

    # The total charge density. 
    # Initiate as np.array to add them
    rho = np.zeros_like(self.z) 

    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
                self.eigenvalues,
                self.eigenstates,
                )
            ):
        # The charge density associated to the nth mode
        # eigenstate is np array. A0 is a callable. 
        # A0(z) is of the same shape as eigenstate
        rho_n = (eigenvalue - self.A0(self.z)) * abs(eigenstate) ** 2
        rho += 1/2*rho_n

    return rho

def calculate_relax_parameter(self):
    """
    In the update law A_{k+1} = (1-c) A_k + c ( - λ (z - 1/2) - ∫ ∫ ρ ), 
    calculate c so that the correction term c( - λ ... ) is exactly self.correction_parameter
    """

    if self.relax_parameter is None:
        c = (
                self.correction_parameter * (- self.lambda_value / 2 + self.a * self.A0_induced(1) )  
                / ( (1-self.correction_parameter) * self.constant_lambda_A0_list[-1](1) 
                    + self.correction_parameter * ( - self.lambda_value / 2 + self.a * self.A0_induced(1))
                    )
                )
    else: 
        c = self.relax_parameter
    
    return c

def calculate_A0_induced(self, rho):
    """
    Calculates the induced A0 as the solution to the Poisson equation.
    It is just integrating the vacuum polarization twice
    """

    rho_interpolated = CubicSpline(self.z, self.rho)
    A0_induced_shifted = rho_interpolated.antiderivative(2)
    offset = A0_induced_shifted(1/2)
    A0_induced = lambda z: -(A0_induced_shifted(z) - offset)

    if A0_induced(1) > 10:
        print(f"n={self.n}, A0_induced(1)={self.A0_induced(1)}")

    return A0_induced

