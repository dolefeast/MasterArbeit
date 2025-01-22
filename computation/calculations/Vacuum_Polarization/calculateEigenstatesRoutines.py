from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from mpmath import quad
import numpy as np


def calculateEigenstates(self):
    """
    Calculates the solutions to the kleinGordon equation by parametrizing the solution family by a parameter, and then looking for the parameter that verifies the boundary conditions
    """

    if self.bcs == "dirichlet":
        # To solve the ODE
        initialValues = (0, 1)
        # To find the zeros when looking for the energies
        # bcsIndex = 0 corresponds to looking for 0 at f(1), 
        # bcsIndex = 1 corresponds to looking for 0 at f'(1), 
        bcsIndex = 0 
    elif self.bcs == "neumann":
        initialValues = (1, 0)
        bcsIndex = 1 
    else:
        NameError(f"{self.bcs} is not a valid boundary condition")

    eigenvalueLowerBound = self.bisectionMethodLowerBound(self.eigenvalues)
    eigenvalueUpperBound = self.bisectionMethodUpperBound(self.eigenvalues)

    parametrizedODE = lambda omega: solve_ivp(
            lambda z, y: self.KleinGordonEquation(z, y, omega), 
            t_span=(0, 1),
            y0=initialValues,
            dense_output=True,
            ) 

    eigenvalues = [
            self.findRootBisection(
                lambda omega: parametrizedODE(omega).sol(1)[bcsIndex], *omegaUpperLower
                ) 
            for omegaUpperLower in zip(
                eigenvalueLowerBound,
                eigenvalueUpperBound
                )
            ]

    # the eigenvalues should be antisymmetric i.e. omegaN = -omega_{-n}
    self.eigenvalues = [ (i-j) / 2 for i, j in zip(eigenvalues, eigenvalues[::-1]) ]
    self.eigenstates = [ parametrizedODE(omega).sol(self.z)[0] for omega in eigenvalues ]

def normalizeEigenstates(self):
    """
    Normalizes the eigenstates 
    """
    # Warning, math

    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
                self.eigenvalues,
                self.eigenstates,
            )
            ):

        eigenstate = CubicSpline(self.z, eigenstate)

        def rhoNWithoutNormalizing(z):
            # Normalizing wrt the symplectic norm
            # the solutions need not be real.
            return (eigenvalue - self.A0(z)) * abs(eigenstate(z))**2

        # Calculate the norm
        normSquared = abs(float(quad(rhoNWithoutNormalizing, [0, 1])))

        norm = self.sign(eigenvalue)*np.sqrt(normSquared)

        self.eigenstates[n] /= norm
