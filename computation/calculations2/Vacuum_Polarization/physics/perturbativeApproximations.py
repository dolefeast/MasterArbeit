from numpy import sin, cos, tan, pi

def perturbativePhi(self, n, z):
    """
    Perturbative mode solution to the TIKGE copied from https://arxiv.org/abs/2010.05499 (equations 16-18)
    """
    omegaN = ( self.m ** 2 + (pi * n) ** 2) ** (1/2)
    if self.bcs == "dirichlet":
        return omegaN ** (-1/2) * (
                sin(pi * n * z) + self.lambdaValue * omegaN / (2 * pi * abs(n)) * (
                        1/(pi * n) * (1/2 - z) * sin(pi*n*z) - z * (1-z) * cos(pi * n * z)
                        )
                    )
    elif self.bcs == "neumann":
        if self.m == 0:
            return ( 2 * self.m ) ** (-1/2) - self.lambdaValue * (2 * m) ** (1/2) * (
                    1/24 - 1/4 * z ** 2 + 1/6 * z ** 3
                    )

        return omegaN ** (-1/2) * (
                cos(pi * n * z) + self.lambdaValue * omegaN / (2 * pi * abs(n)) * (
                        1/(pi * n) * (1/2 - z) * cos(pi*n*z) + (z * (1-z) + (pi*n)**(-2)) * sin(pi * n * z)
                        )
                    )

def perturbativeModeRho(self, n, z):
    """
    Perturbative charge density of the mode N copied from https://arxiv.org/abs/2010.05499 
    """
    
    assert self.bcs == "dirichlet", "Only the Dirichlet boundary conditions have been contemplated"

    return - 2 * self.lambdaValue * (
            (z - 1/2) * (abs(omegaN)/(pi*n)**2 - 1/abs(omegaN)) * sin(n * pi * z)**2
            + abs(omegaN) / pi / n * sin(pi * n * z) * cos(pi * n * z) * (1-z) * z
            )

def cotZ(z):
    """
    The z*(1-z)cot(pi * z) appearing in equation 38 of https://arxiv.org/abs/2010.05499.
    The z = 0, 1 exist in the limit, which are calculated using L'Hôpital's rule
    """

    if z == 0:
        return -1/pi
    elif z == 1:
        return 1/pi
    else:
        return z * (z - 1) * tan(pi * z) ** -1

def perturbativeTotalVacuumPolarization(self, n, z):
    """
    Perturbative vacuum polarization of the Klein-Gordon field after doing the summation (for the massless case). Copied from https://arxiv.org/abs/2010.05499 
    """
    assert self.bcs == "dirichlet", "Only the Dirichlet boundary conditions have been contemplated"

    return - self.e * self.lambdaValue * cotZ(z) / 2
