from numpy import shape
from scipy.interpolate import CubicSpline

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

    # Theres a minus missing!
    A0_induced = E_induced.antiderivative(1)
    
    # It is added here
    self.A0_induced = -A0_induced(self.z) + A0_induced(1/2)
