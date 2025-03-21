import numpy as np

def KleinGordonEquation(self, z, y, omega): 
    """
    The time-independent Klein-Gordon equation. This is input into the ODE solver 
    """
    backgroundField = self.A0(z)

    kleinGordonArray = np.array(
            (y[1],
                - (( omega - backgroundField) ** 2 - self.m ** 2) * y[0] )
            )

    return kleinGordonArray
