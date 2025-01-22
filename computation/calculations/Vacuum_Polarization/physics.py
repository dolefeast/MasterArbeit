import numpy as np

def eigenvalueGuess(self, n):
    """The λ=0 eigenvalue generator for mode n """
    return self.sign(n) * np.sqrt( (n * np.pi)**2 + self.m ** 2 )

def KleinGordonEquation( self, z, y, omega):
    """
    The equation of motion of the field
    """
    # The differential equation

    backgroundField = self.A0(z)

    kleinGordon = np.array(
        (y[1], -((omega - backgroundField) ** 2 + self.m ** 2) * y[0])
    )

    return kleinGordon
