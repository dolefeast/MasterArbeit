import numpy as np

def perturbativeEigenvalue(self, n):
    """The λ=0 eigenvalue generator for mode n """
    return self.sign(n) * np.sqrt( (n * np.pi)**2 + self.m ** 2 )

def freePhi(self, n):
    omega = self.perturbativeEigenvalue(n)
    if self.bcs == "neumann":
        if n==0:
            return lambda z: (2 * self.m) ** -1/2
        else:
            return lambda z: abs(omega) ** -1/2 * np.cos(n*np.pi*z)
    elif self.bcs == "dirichlet":
        return lambda z: abs(omega) ** -1/2 * np.sin(n*np.pi*z)

def perturbativePhi(self, n:int) -> callable:
    """
    Returns a function of z corresponding to the mode N, lambda value lambdaValue and mass m
    """
    if self.bcs == "neumann":
        if n == 0:
            return lambda z: (2*self.m) ** -1/2 - self.lambdaValue * np.sqrt(2*self.m) * (
                    1/24 - 1/4 * z ** 2 + 1/6 * z ** 3 
                    )
        omega = self.perturbativeEigenvalue(n)
        return lambda z: abs(omega) ** -1/2 * (
                np.cos(np.pi * n * z)
                + self.lambdaValue * omega / 2 / np.pi / n * 
                (1 / np.pi / n * (1/2-z) * np.cos(np.pi * n * z) + 
                    (
                        z * (1-z) + (np.pi * n)**-2) *np.sin(np.pi * n * z)
                        )
                )
    elif self.bcs == "dirichlet":
        if n == 0:
            raise ValueError("Dirichlet boundary conditions does not admit n=0")
        omega = self.perturbativeEigenvalue(n)
        return lambda z: abs(omega) ** -1/2 * (
                np.sin(np.pi * n * z)
                + self.lambdaValue * omega / 2 / np.pi / n * 
                (1 / np.pi / n * (1/2-z) * np.sin(np.pi * n * z) - 
                    (
                        z * (1-z) ) *np.cos(np.pi * n * z)
                        )
                )

def KleinGordonEquation( self, z, y, omega):
    """
    The equation of motion of the field
    """
    # The differential equation

    backgroundField = self.A0(z)

    kleinGordon = np.array(
        (y[1], -((omega - backgroundField) ** 2 - self.m ** 2) * y[0])
    )

    return kleinGordon

