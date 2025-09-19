import matplotlib.pyplot as plt
import numpy as np

from Vacuum_Polarization import Vacuum_Polarization

fig, axPhi = plt.subplots()
fig, axRho = plt.subplots()

vp = Vacuum_Polarization(maxN = 10, lambdaMin=1., bisectionTol=1e-15)

# vp.calculateEigenstatesParallel()
vp.eigenstates = [
vp.calculateSingleEigenstate(vp.eigenvalues[-3:-1])[1] ]

n = vp.maxN + -1
omegaN = vp.eigenvalues[-1]
phiN = np.array(vp.eigenstates[-1])
axPhi.plot(vp.z, phiN, label="non perturbative")
axPhi.plot(vp.z, vp.perturbativePhi(n, vp.z), label="perturbative")

vp.A0 = lambda z : -vp.lambdaValue * (z - 1/2)

rho = (
        (omegaN - vp.A0(vp.z)) * phiN**2 
        + (- omegaN - vp.A0(vp.z)) * phiN[::-1]**2
        )
# vp.calculateRho()
# rho = vp.rho

axRho.plot(vp.z, rho, label="non perturbative")
axRho.plot(vp.z, vp.perturbativeModeRho(n, vp.z), label="perturbative")

axPhi.set_ylabel(f"$\phi_{n}$")
axPhi.set_xlabel(f"z")
axPhi.legend()

axRho.set_ylabel(r"$\rho_{n}$")
axRho.set_xlabel(f"z")
axRho.legend()

plt.show()
