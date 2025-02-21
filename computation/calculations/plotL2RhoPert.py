import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import numpy as np
import scienceplots

plt.style.use(["science", "high-contrast"])

from plotScripts.readDataScripts import openSolutionFamilyArray
from Vacuum_Polarization import Vacuum_Polarization

def maxDerivative(y):
    x = np.linspace(0, 1, len(y))
    y = CubicSpline(x, y)

    return -y.derivative()(1/2)


directory, m, a, solutionFamilyArray = openSolutionFamilyArray("")

fig1, axMaxRho = plt.subplots()
fig2, axCharge = plt.subplots()
fig3, axChargeComparison = plt.subplots()

maxRhoList = [max(rho) for rho in solutionFamilyArray["rho"]]
maxPertRhoList = []
for lambdaValue in solutionFamilyArray["lambdaValue"]:
    vp = Vacuum_Polarization(m=0, a=1, lambdaMin=lambdaValue)

    maxPertRhoList.append(max(vp.perturbativeVacuumPolarizationMasslessInf(vp.z) - vp.e**2 / np.pi * lambdaValue * (vp.z - 1/2)))
vp = Vacuum_Polarization(m=0, a=0, lambdaMin=1)

maxQList = [maxDerivative(A0) for A0 in solutionFamilyArray["A0Induced"]]

axMaxRho.plot(solutionFamilyArray["lambdaValue"], maxRhoList, label = r"Non perturbative max($\rho$)")
axMaxRho.plot(solutionFamilyArray["lambdaValue"], maxPertRhoList, label = r"Perturbative max($\rho$)")
axCharge.plot(solutionFamilyArray["lambdaValue"], maxQList)
axChargeComparison.plot(solutionFamilyArray["lambdaValue"], np.array(solutionFamilyArray["lambdaValue"]) + np.array(maxQList))

axMaxRho.set_xlabel(r"$\lambda$")
axCharge.set_xlabel(r"$\lambda$")
axChargeComparison.set_xlabel(r"$\lambda$")

axMaxRho.set_ylabel(r"max($\rho$)")
axCharge.set_ylabel(r"$Q_I$")
axChargeComparison.set_ylabel(r"$\lambda + Q_I$")


fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

plt.show()
