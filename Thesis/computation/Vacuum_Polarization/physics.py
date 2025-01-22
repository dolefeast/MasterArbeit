import numpy as np

def eigenvalue_guess(self, n):
    """The λ=0 eigenvalue generator for mode n """
    return self.sign(n) * np.sqrt( (n * np.pi)**2 + self.m ** 2 )

def Klein_Gordon_equation( self, z, y, omega):
    """
    The equation of motion of the field
    """
    # The differential equation

    background_field = self.A0(z)

    klein_gordon = np.array(
        (y[1], -((omega - background_field) ** 2 + self.m ** 2) * y[0])
    )

    return klein_gordon
