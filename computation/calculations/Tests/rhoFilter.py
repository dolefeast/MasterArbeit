from Vacuum_Polarization import Vacuum_Polarization
from plotScripts.readDataScripts import openSolutionFamilyArray

import matplotlib.pyplot as plt
import numpy as np

fig1, ax1 = plt.subplots(figsize=(10, 10))

def chargeDensityComparisonPlot(compute):

    fig2, ax2 = plt.subplots(figsize=(10, 10))
    fig2.suptitle(f"$\lambda = {compute.lambdaValue}$, comparing non perturbative and perturbative" + r" $\rho_1$")
    ax2.plot((compute.eigenvalues[compute.maxN] - compute.e * compute.A0(compute.z)) * compute.eigenstates[compute.maxN]**2 + 0.01, label=r"Non perturbative $\rho$")
    ax2.plot((compute.perturbativeEigenvalue(1) + compute.lambdaValue * (compute.z - 1/2)) * compute.perturbativePhi(1)(compute.z)**2, label=r"perturbative $\rho$")

    ax2.legend(loc="best")

def filteringComparisonPlot(compute, ax1):
    # This should be the unfiltered vacuum polarization
    # One adds the potential at the end because the sum of the mode vacuum polarizations
    # tends to "rest" over the curve e²/pi A0(z)
    unfilteredRhoSubtraction = compute.calculateRho() - compute.e**2 / np.pi * compute.lambdaValue * (compute.z -  1/2) 
    unfilteredRhoSubtraction = antisymmetrize(unfilteredRhoSubtraction)

    x, filteredRho = compute.extendAndFilter(compute.rho)
    filteredRho = antisymmetrize(filteredRho)
    x, filteredRhoSubtraction = compute.extendAndFilter(unfilteredRhoSubtraction)

    # ax1.plot(unfilteredRhoSubtraction - compute.rho)
        # , label= r"$\sum_{n=0}^{"+str(compute.maxN)+r"} \rho_n^m - \rho_n^0  -  \frac{e\lambda}{2} z(1-z) \cot(\pi z)$")


    # ax1.plot(filteredRho, label=r"filtered $\rho$ without subtraction an addition of massless terms")
    ax1.plot(compute.z, filteredRhoSubtraction, label=r"Filtered $\rho$ with subtraction, N = " + str(compute.maxN))
    #, label=
    #       #r"filtered $-\sum\rho^0 - \frac{e\lambda}{2} z(z-1) \cot \pi z$ ")
    #       r"filtered $\sum(\rho^{True}_n - \rho^{pert, massless}_n) - \frac{e\lambda}{2} z(z-1) \cot\pi z$")


def antisymmetrize(arr:list):
    return [(xForward - xBackwards)/2 for xForward, xBackwards in zip(arr, arr[::-1])]

for _ in range(4):
    filterRegex = ""

    directory, solutionFamilyArray = openSolutionFamilyArray(filterRegex)

    solutionFamily = {key:value[-1] for key, value in solutionFamilyArray.items()}

    compute = Vacuum_Polarization(bcs=directory.name, subtractMasslessPertVacuumPolarization=True)
    compute.setConfigFromDict(solutionFamily)

    filteringComparisonPlot(compute, ax1)

ax1.legend(loc="best")

fig1.suptitle(f"$\lambda = {compute.lambdaValue}$, with Wernerssons+Zahn's convolution")

ax1.set_xlabel("z")
plt.show()
