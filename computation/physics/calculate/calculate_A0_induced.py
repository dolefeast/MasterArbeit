from numpy import shape
from scipy.interpolate import CubicSpline
from physics.utils.Fields import Field

def calculate_A0_induced(self):
    """
    Calculates the induced A0_field from the current self.rho
    """
    assert self.rho.ndim == 1, "self.rho.dim should be 1-dimensional array"
    rho = CubicSpline(
            self.z,
            self.rho
            )

    E_induced = rho.antiderivative(1)
    self._E_induced = E_induced(self.z)

    # There is a minus missing!
    A0_induced = E_induced.antiderivative(1)
    
    # It is added here
    self.A0_induced = -A0_induced(self.z) + A0_induced(1/2)
    self.A0_field = Field(
            n_points = self.n_points,
            value = -self.lambda_value * (self.z - 1/2) + self.a * self.A0_induced,
            )
