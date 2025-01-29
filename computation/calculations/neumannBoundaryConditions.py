from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

import scienceplots
plt.style.use(["science"])

from Vacuum_Polarization import Vacuum_Polarization

def KleinGordonEquation(z, y, omega, lambdaValue):
    A0 = -lambdaValue*(z-1/2)

    kleinGordon = np.array(
        (y[1], -((omega - backgroundField) ** 2 + self.m ** 2) * y[0])
    )

def firstModeOmegaAsFunctionOfMass(massMin:float, massMax:float, massNPoints:int, lambdaValue:float, bcs:str)-> list:
    massRange = np.linspace(massMin, massMax, massNPoints)
    omegaRange = [main(m, lambdaValue, bcs)[0] for m in massRange]

    return massRange, omegaRange

def plotOmegaVsMass(*args)->None:
    plt.plot(*firstModeOmegaAsFunctionOfMass(*args)
            )

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

def main(m:float, lambdaValue:float, bcs:str)->float:
    compute = Vacuum_Polarization(m=m, lambdaMin=lambdaValue, bcs=bcs)

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
        omega = compute.findRootBisection(lambda omega: parameterizedODE(omega).sol(1)[bcsIndex], 0, m+0.1)
    except RuntimeError:
        return None, None

    return omega, parameterizedODE



fig, ax = plt.subplots(figsize=(10,8))

z = np.linspace(0, 1, 200)
m = 6
lambdaValue = 0.1
bcs = "neumann"

omega, parameterizedODE = main(m, lambdaValue, bcs)
secondMode = perturbativePhi(0, lambdaValue, m, bcs=bcs)

plt.plot(z, parameterizedODE(omega).sol(z)[0]/max(parameterizedODE(omega).sol(z)[0]), label="Numerical solution")
plt.plot(z, secondMode(z)/max(secondMode(z)), label="Perturbative solution")

plt.legend(loc="best")
# plt.suptitle(f"m = {m}, $\lambda = {lambdaValue}$, $\omega = {round(omega, 2)}$")

# print( plotOmegaVsMass(0, 10, 50, 0.1, "dirichlet"))

ax.set_xlabel("z")
plt.show()
