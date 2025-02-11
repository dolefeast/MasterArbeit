import numpy as np

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

