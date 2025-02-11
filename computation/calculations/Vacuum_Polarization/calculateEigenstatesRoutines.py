from scipy.integrate import solve_ivp
from functools import partial
# from multiprocessing import Pool
# from concurrent.futures import ProcessPoolExecutor


from scipy.interpolate import CubicSpline
from mpmath import quad
import numpy as np


def calculateNormSquared(self, omega, eigenstate):
    normSquared = abs(float(
        quad(lambda z: 2 * ( omega - self.e * self.A0(z) ) * eigenstate.sol(z)[0] ** 2, [0, 1])
        ))

    return normSquared

def calculateSingleEigenstate(self, omegaUpperLower):
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

    # This IVP solution propagates the initial conditions at the left boundary to the right
    # without knowing where it is going to land
    parameterizedODE = lambda omega: solve_ivp(
            lambda z, y: self.KleinGordonEquation(z, y, omega), 
            t_span=(0, 1),
            y0=initialValues,
            dense_output=True,
            ) 

    # bcsIndex indicates with respect to which derivative of the field is the omega found
    # If bcs=dirichlet we want the zeroth derivative to be zero at the right boundary (z=1)
    # If bcs=neumann we want the first derivative to be zero at the right boundary (z=1)
    omega = self.findRootBisection(lambda omega:parameterizedODE(omega).sol(1)[bcsIndex], *omegaUpperLower)

    normCalc = self.calculateNormSquared(omega, parameterizedODE(omega))
    return omega, parameterizedODE(omega).sol(self.z)[0]/np.sqrt(abs(normCalc))

def calculateEigenstatesParallel(self):
    from pathos.multiprocessing import ProcessingPool as Pool

    p = Pool(4)

    def _calculateSingleEigenstate(omegaUpperLower):
        return calculateSingleEigenstate(self, omegaUpperLower)

    omegaUpperArray = self.bisectionMethodUpperBound(self.eigenvalues)
    omegaLowerArray = self.bisectionMethodLowerBound(self.eigenvalues)

    eigenValuesEigenStates = p.map(_calculateSingleEigenstate, zip(omegaUpperArray, omegaLowerArray))
    # with Pool() as pool:
    #     eigenValuesEigenStates = pool.map(_calculateSingleEigenstate, zip(omegaLowerArray, omegaUpperArray))

    self.eigenvalues = [ element[0] for element in eigenValuesEigenStates ]
    self.eigenstates = [ element[1] for element in eigenValuesEigenStates ]

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

    parameterizedODE = lambda omega: solve_ivp(
            lambda z, y: self.KleinGordonEquation(z, y, omega), 
            t_span=(0, 1),
            y0=initialValues,
            dense_output=True,
            ) 

    eigenvalues = [
            self.findRootBisection(
                lambda omega: parameterizedODE(omega).sol(1)[bcsIndex], *omegaUpperLower
                ) 
            for omegaUpperLower in zip(
                eigenvalueLowerBound,
                eigenvalueUpperBound
                )
            ]

    # the eigenvalues should be antisymmetric i.e. omegaN = -omega_{-n}
    self.eigenvalues = [ (i-j) / 2 for i, j in zip(eigenvalues, eigenvalues[::-1]) ]
    self.eigenstates = [ parameterizedODE(omega).sol(self.z)[0] for omega in eigenvalues ]

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
            return 2 * (eigenvalue - self.A0(z)) * abs(eigenstate(z))**2

        # Calculate the norm
        normSquared = abs(float(quad(rhoNWithoutNormalizing, [0, 1])))

        norm = self.sign(eigenvalue)*np.sqrt(normSquared)

        self.eigenstates[n] /= norm
