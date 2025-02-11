import numpy as np

def perturbativeEigenvalue(self, n):
    """The λ=0 eigenvalue generator for mode n """
    return self.sign(n) * np.sqrt( (n * np.pi)**2 + self.m ** 2 )

def freePhi(self, n):
    omega = self.perturbativeEigenvalue(n)
    if self.bcs == "neumann":
        if n==0:
            return lambda z: (2 * self.m) ** -.5
        else:
            return lambda z: abs(omega) ** -.5 * np.cos(n*np.pi*z)
    elif self.bcs == "dirichlet":
        return lambda z: abs(omega) ** -.5 * np.sin(n*np.pi*z)

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
        return lambda z: abs(omega) ** -0.5 * (
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
        return lambda z: abs(omega) ** -0.5 * (
                np.sin(np.pi * n * z)
                + self.lambdaValue * omega / 2 / np.pi / n * 
                (1 / np.pi / n * (1/2-z) * np.sin(np.pi * n * z) - 
                    (
                        z * (1-z) ) *np.cos(np.pi * n * z)
                        )
                )

def sign(n):
    if n > 0:
        return 1
    elif n < 0:
        return -1
    return 0

def perturbativeVacuumPolarizationMasslessInf(self, z) -> float:
    # The evaluation of z (z-1) * cot(piz) at z=0,1 is not defined, but the limit is (L'Hôpital)
    limitAtZ0 = -1/np.pi
    limitAtZ1 = 1/np.pi
    zReduced = z[1:-1]
    cotArray = np.append(limitAtZ0, np.tan(np.pi* zReduced)**-1*zReduced*(1-zReduced))
    cotArray = np.append(cotArray, limitAtZ1)
    return -1/2* self.e * self.lambdaValue * cotArray

def perturbativeVacuumPolarizationNModeMassless(self, n, z) -> float:
    return -  self.e * self.lambdaValue *  np.sin(2 * n*np.pi*z)  * (1-z) * z

