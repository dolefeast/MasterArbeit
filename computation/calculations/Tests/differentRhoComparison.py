from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
from mpmath import quad

import scienceplots
# plt.style.use(["science"])

from Vacuum_Polarization import Vacuum_Polarization

def firstModeOmegaAsFunctionOfMass(massMin:float, massMax:float, massNPoints:int, lambdaValue:float, bcs:str)-> list:
    massRange = np.linspace(massMin, massMax, massNPoints)
    omegaRange = [calculateEigens(m, lambdaValue, bcs)[0] for m in massRange]

    return massRange, omegaRange

def plotOmegaVsMass(*args)->None:
    plt.plot(*firstModeOmegaAsFunctionOfMass(*args)
            )

def calculateNorm(omega, eigenstate, lambdaValue):

    try:
        normSquared = abs(float(quad(lambda z: 2* (omega + lambdaValue * (z - 1/2)) * eigenstate.sol(z)[0]**2, [0, 1] )))
        return normSquared
    except AttributeError:
        eigenstate = CubicSpline(compute.z, eigenstate)

        normSquared = abs(float(quad(lambda z: 2* (omega + lambdaValue * (z - 1/2)) * eigenstate(z)**2, [0, 1] )))
        return normSquared

def calculateEigens(compute:object, n:int, window=np.pi/2)->float:

    if compute.bcs == "neumann":
        initialValues = (1, 0)
        bcsIndex = 1
    elif compute.bcs == "dirichlet":
        initialValues = (0, 1)
        bcsIndex = 0

    parameterizedODE = lambda omega: solve_ivp(lambda z, y: compute.KleinGordonEquation(z, y, omega),
            t_span=(0, 1),
            y0=initialValues, 
            dense_output=True)

    try:
        omega = compute.findRootBisection(lambda omega: parameterizedODE(omega).sol(1)[bcsIndex], compute.perturbativeEigenvalue(n)-window/2, compute.perturbativeEigenvalue(n)+window/2)
    except RuntimeError:
        return None, None
    
    normCalc = calculateNorm(omega, parameterizedODE(omega), compute.lambdaValue)

    return omega, parameterizedODE(omega).sol(compute.z)[0]/np.sqrt(abs(normCalc))

def vacuumPolarizationComparison(fig1, ax1, compute, n):
    # Non perturbative charge densities

    assert n>0, "n must be a positive integer"

    compute.calculateEigenstates()
    compute.normalizeEigenstates()

    nonPertRhoPositive =  (compute.eigenvalues[compute.maxN + n - 1] + compute.lambdaValue * (compute.z-1/2)) * (compute.eigenstates[compute.maxN + n - 1])**2
    nonPertRhoNegative =  (compute.eigenvalues[compute.maxN-n] + compute.lambdaValue * (compute.z-1/2)) * (compute.eigenstates[compute.maxN-n])**2

    totalnonPertRho = nonPertRhoPositive + nonPertRhoNegative

    # Perturbative vacuum polarization
    pertRho = lambda z: (omega + lambdaValue * (z-1/2)) * (perturbativeMode(z))**2
    oppositePertRho = lambda z: -pertRho(1-z)
    totalPertRho = lambda z: pertRho(z) + oppositePertRho(z)

    # Free field vacuum polarization
    freeRho = lambda z: (omega + lambdaValue * (z-1/2) ) * (freePhi(z))**2
    oppositeFreeRho = lambda z: -freeRho(1-z)
    totalFreeRho = lambda z: freeRho(z) + oppositeFreeRho(z)

    # ax1.plot(z, totalnonPertRho(z)/max(totalnonPertRho(z)), label=r"Numerical $\rho"+"_{"+str(n)+"}$")
    ax1.plot(compute.z, totalnonPertRho, label=r"Numerical $\rho"+"_{"+str(n)+"}$")
    ax1.plot(compute.z, totalPertRho(compute.z), label=r"Perturbative $\rho"+"_{"+str(n)+"}$")
    ax1.plot(compute.z, totalFreeRho(compute.z), label=r"Free $\rho"+"_{"+str(n)+"}$")
    # ax1.plot(z, freePhi(z)/max1(freePhi(z)), label="Free solution")

    fig.suptitle("$\phi$ comparison")

def fieldSolutionComparisons(fig, ax, omega:float, parameterizedODE:callable):
    ax.plot(z, parameterizedODE, label="Numerical solution")
    ax.plot(z, perturbativeMode(z), label="Perturbative solution")
    ax.plot(z, freePhi(z), label="Free solution")

    fig.suptitle("$\phi$ comparison")


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(10,8))
    fig1, ax1 = plt.subplots(figsize=(10,8))


    n=1
    z = np.linspace(0, 1, 100*(n+1))
    m = 0
    lambdaValue = 1

    bcs = "dirichlet"

    compute = Vacuum_Polarization(m=m, lambdaMin=lambdaValue, bcs=bcs, maxN = n)
    omega, parameterizedODE = calculateEigens(compute, n )
    perturbativeMode = compute.perturbativePhi(n)
    freePhi = compute.freePhi(n)

    fieldSolutionComparisons(fig, ax, omega, parameterizedODE)
    vacuumPolarizationComparison(fig1, ax1, compute, n)


    # ax1.plot(z, parameterizedODE(omega).sol(z)[0]/max(parameterizedODE(omega).sol(z)[0]), label="Numerical solution")
    # ax1.plot(z, perturbativeMode(z)/max(perturbativeMode(z)), label="Perturbative solution")
    # ax1.plot(z, freePhi(z)/max(freePhi(z)), label="Free solution")

    ax.legend(loc="best")
    ax1.legend(loc="best")
    # plt.suptitle(f"m = {m}, $\lambda = {lambdaValue}$, $\omega = {round(omega, 2)}$")

    # print( plotOmegaVsMass(0, 10, 50, 0.1, "dirichlet"))

    ax.set_xlabel("z")
    ax1.set_xlabel("z")

    fig.suptitle(r"$\phi_{"+str(n)+"}$ comparison, $\lambda=$" + str(compute.lambdaValue) )
    fig1.suptitle(r"$\rho_{"+str(-n)+r"} + \rho_{"+str(n)+"}$ comparison, $\lambda=$" + str(compute.lambdaValue) )

    plt.show()
