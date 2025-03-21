from scipy.integrate import solve_ivp
from scipy.integrate import quad as spquad
from mpmath import quad, quadsubdiv
import mpmath
import numpy as np

mpmath.mp.dps = 40
mpmath.mp.prec = 50

def calculateNormSquared(self, omega: float, eigenstate: callable, error=False, verbose=False):
    """
    Calculates the squared norm of a calculated non-normalized eigenstate
    """

    normSquared = quad(lambda z: 2 * ( omega - self.e * self.A0(z) ) * eigenstate.sol(z)[0] ** 2, [0, 1], error=error, verbose=verbose) 

    if error:
        print(f"mpmath.quad calculates norm = {normSquared[0]} ± {100*normSquared[1]/2/normSquared[0]}%")

    return normSquared

def calculateSingleEigenstate(self, eigenvalueRange: tuple, error=False, verbose=False):
    """
    Calculates a single eigenstate of the KG field with eigenvalue in the range (eigenvalueRange[0], eigenvalueRange[1])
    """

    if self.bcs == "dirichlet":
        # Initial values of the IVP
        initialValues = (0, 1)
        # To find the roots of phi[bcsIndex] when applying the boundary conditions
        bcsIndex = 0 
    elif self.bcs == "neumann":
        initialValues = (1, 0)
        bcsIndex = 1 

    parameterizedODE = lambda omega: solve_ivp(
            lambda z, y: self.KleinGordonEquation(z, y, omega), 
            t_span=(0, 1),
            y0=initialValues,
            dense_output=True,
            rtol=1e-10,
            atol=1e-10,
            )

    eigenvalue = self.bisectionMethod(
            lambda omega: parameterizedODE(omega).sol(1)[bcsIndex],
            *eigenvalueRange,
            maxiter=self.bisectionMaxiter,
            tol=self.bisectionTol,
            )

    normSquared = self.calculateNormSquared(eigenvalue, parameterizedODE(eigenvalue), error=error, verbose=verbose)

    return eigenvalue, parameterizedODE(eigenvalue).sol(self.z)[0]/np.sqrt(abs(normSquared))

def calculateEigenstatesParallel(self):
    if self.parallelization:
        from pathos.multiprocessing import ProcessingPool as Pool

        p = Pool(12)

    def _calculateSingleEigenstate(omegaUpperLower):
        return calculateSingleEigenstate(self, omegaUpperLower)

    omegaUpperArray = self.eigenvaluesUpperBound(self.eigenvalues)
    omegaLowerArray = self.eigenvaluesLowerBound(self.eigenvalues)

    if self.parallelization:
        eigenValuesEigenStates = p.map(_calculateSingleEigenstate, zip(omegaUpperArray, omegaLowerArray))
    else: 
        raise NotImplementedError("Parallelization must stay on")
        eigenValuesEigenStates = map(_calculateSingleEigenstate, zip(omegaUpperArray, omegaLowerArray))


    self.eigenvalues = [ element[0] for element in list(eigenValuesEigenStates) ]
    self.eigenstates = [ element[1] for element in list(eigenValuesEigenStates) ]
