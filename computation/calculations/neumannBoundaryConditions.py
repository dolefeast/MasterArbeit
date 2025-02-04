from scipy.integrate import solve_ivp
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

def calculateNorm(eigenvalue, eigenstate, lambdaValue):
    try:
        normSquared = abs(float(quad(lambda z: (eigenvalue + lambdaValue * (z - 1/2)) * eigenstate.sol(z)[0]**2, [0, 1] ))) 
    except AttributeError:
        normSquared = abs(float(quad(lambda z: (eigenvalue + lambdaValue * (z - 1/2)) * eigenstate(z)**2, [0, 1] ))) 

    return normSquared

def perturbativePhi(n:int, lambdaValue:float, m:float, bcs:str) -> callable:
    """
    Returns a function of z corresponding to the mode N, lambda value lambdaValue and mass m
    """
    if bcs == "neumann":
        if n == 0:
            return lambda z: (2*m) ** -1/2 - lambdaValue * np.sqrt(2*m) * (
                    1/24 - 1/4 * z ** 2 + 1/6 * z ** 3 
                    )
        omega = abs(n)/n * (m**2 + np.pi**2 * n**2)**1/2
        return lambda z: abs(omega) ** -1/2 * (
                np.cos(np.pi * n * z)
                + lambdaValue * omega / 2 / np.pi / n * 
                (1 / np.pi / n * (1/2-z) * np.cos(np.pi * n * z) + 
                    (
                        z * (1-z) + (np.pi * n)**-2) *np.sin(np.pi * n * z)
                        )
                )
    elif bcs == "dirichlet":
        if n == 0:
            raise ValueError("Dirichlet boundary conditions does not admit n=0")
        omega = abs(n)/n * (m**2 + np.pi**2 * n**2)**1/2
        return lambda z: abs(omega) ** -1/2 * (
                np.sin(np.pi * n * z)
                + lambdaValue * omega / 2 / np.pi / n * 
                (1 / np.pi / n * (1/2-z) * np.sin(np.pi * n * z) - 
                    (
                        z * (1-z) ) *np.cos(np.pi * n * z)
                        )
                )

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

    return omega, parameterizedODE, normCalc

def vacuumPolarizationComparison(fig1, ax1, compute, n):
    # Non perturbative vacuum polarizations
    nonPertOmegaPositive, parameterizedODEPositive, norm = calculateEigens(compute, n)
    nonPertRhoPositive = lambda z: (nonPertOmegaPositive + compute.lambdaValue * (z-1/2)) * (parameterizedODEPositive(nonPertOmegaPositive).sol(z)[0])**2/norm

    nonPertOmegaNegative, parameterizedODENegative, norm = calculateEigens(compute, -n)
    nonPertRhoNegative = lambda z: (nonPertOmegaNegative + compute.lambdaValue * (z-1/2)) * (parameterizedODENegative(nonPertOmegaNegative).sol(z)[0])**2/norm

    totalnonPertRho = lambda z: nonPertRhoPositive(z) + nonPertRhoNegative(z)

    # Perturbative vacuum polarizations
    pertRho = lambda z: (omega + lambdaValue * (z-1/2)) * (perturbativeMode(z))**2
    oppositePertRho = lambda z: -pertRho(1-z)
    totalPertRho = lambda z: pertRho(z) + oppositePertRho(z)

    # Free field vacuum polarizations
    freeRho = lambda z: (omega ) * (freePhi(z))**2
    oppositeFreeRho = lambda z: -freeRho(1-z)
    totalFreeRho = lambda z: freeRho(z) + oppositeFreeRho(z)

    # ax1.plot(z, totalnonPertRho(z)/max(totalnonPertRho(z)), label=r"Numerical $\rho"+"_{"+str(n)+"}$")
    # ax1.plot(z, totalPertRho(z)/max(totalPertRho(z)), label=r"Perturbative $\rho"+"_{"+str(n)+"}$")
    # ax1.plot(z, totalFreeRho(z)/max(totalFreeRho(z)), label=r"Free $\rho"+"_{"+str(n)+"}$")
    ax1.plot(z, totalnonPertRho(z)/max(totalnonPertRho(z)), label=r"Numerical $\rho"+"_{"+str(n)+"}$")
    ax1.plot(z, totalPertRho(z)/max(totalPertRho(z)), label=r"Perturbative $\rho"+"_{"+str(n)+"}$")
    ax1.plot(z, totalFreeRho(z), label=r"Free $\rho"+"_{"+str(n)+"}$")
    # ax1.plot(z, freePhi(z)/max1(freePhi(z)), label="Free solution")

    fig.suptitle("$\phi$ comparison")
def fieldSolutionComparisons(fig, ax, omega:float, parameterizedOde:callable):
    ax.plot(z, parameterizedODE(omega).sol(z)[0]/max(parameterizedODE(omega).sol(z)[0]), label="Numerical solution")
    ax.plot(z, perturbativeMode(z)/max(perturbativeMode(z)), label="Perturbative solution")
    ax.plot(z, freePhi(z)/max(freePhi(z)), label="Free solution")

    fig.suptitle("$\phi$ comparison")


fig, ax = plt.subplots(figsize=(10,8))
fig1, ax1 = plt.subplots(figsize=(10,8))

n=50
z = np.linspace(0, 1, 100*(n+1))
m = 6
lambdaValue = 4

bcs = "dirichlet"

compute = Vacuum_Polarization(m=m, lambdaMin=lambdaValue, bcs=bcs)

omega, parameterizedODE, norm = calculateEigens(compute, n )
perturbativeMode = compute.perturbativePhi(n)
freePhi = compute.freePhi(n)

fieldSolutionComparisons(fig, ax, omega, parameterizedODE)
vacuumPolarizationComparison(fig1, ax1, compute, n)

# ax1.plot(z, parameterizedODE(omega).sol(z)[0]/max(parameterizedODE(omega).sol(z)[0]), label="Numerical solution")
# ax1.plot(z, perturbativeMode(z)/max(perturbativeMode(z)), label="Perturbative solution")
# ax1.plot(z, freePhi(z)/max(freePhi(z)), label="Free solution")

plt.legend(loc="best")
# plt.suptitle(f"m = {m}, $\lambda = {lambdaValue}$, $\omega = {round(omega, 2)}$")

# print( plotOmegaVsMass(0, 10, 50, 0.1, "dirichlet"))

ax.set_xlabel("z")
plt.show()
