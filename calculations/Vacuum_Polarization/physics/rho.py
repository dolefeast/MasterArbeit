import numpy as np
from mpmath import mpf


def calculateRho(self):
    """
    The vacuum polarization corresponding to the current family of eigenstates 
    """
    rho = np.array([mpf("0.")] * self.nPoints)

    for n, (omegaN, phiN) in enumerate(
            zip(
                self.eigenvalues,
                self.eigenstates,
                )
            ):
        # The charge density associated to the nth mode
        # eigenstate is np array. A0 is a callable. 
        # A0(z) is of the same shape as eigenstate
        rhoN = self.a*self.e*((omegaN - self.e*self.a*self.A0(self.z)) * phiN ** 2 
                + (-omegaN - self.e*self.a*self.A0(self.z)) * phiN[::-1] ** 2)


        if self.subtractPertModes:
            start = self.bcs=="dirichlet"
            rhoN -= self.perturbativeModeRho(n+start, self.z)
        rho += np.array(rhoN)

    if self.subtractPertModes:
        self.rho = rho + self.perturbativeTotalVacuumPolarization(self.z)
    else:
        self.rho = rho

    if self.hadamard:
        print("Adding extra potential term")
        self.rho += self.e**2 * self.a**2 / np.pi * self.A0(self.z)

def convolveRho(self, kernel = None):
    from scipy.interpolate import CubicSpline

    if kernel is None:
        window = self.nPoints // (self.maxN + 1) *2
        kernel = [ np.exp(-x**2*4/window**2) for x in range(-3*window , 3*window+1)]

    # Always normalize the convolution kernel
    convolutionSum = sum(kernel)
    kernel = [i/convolutionSum for i in kernel]

    convolvedRho = np.convolve(self.rho, kernel, mode="same")

    if self.bcs == "dirichlet" and True:
    # The indices outside of the boundary effects
        idx = np.where(np.logical_and(
            self.z>window/self.nPoints,
            self.z<1-window/self.nPoints)
            )

        zReduced = self.z[idx]
        rhoReduced = convolvedRho[idx]

        zReduced = np.concatenate(([0], zReduced, [1]))
        rhoReduced = np.concatenate(([0], rhoReduced, [0]))

        convolvedRho = CubicSpline(zReduced, rhoReduced)(self.z)
    
    self.rho = convolvedRho
