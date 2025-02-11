from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
from mpmath import quad

from Vacuum_Polarization import Vacuum_Polarization
from plotScripts.readFiles import getDirectoryMA
from plotScripts.readDataScripts import getPosixForQuantities, openPosixDict, strToFloat

fig1, ax1 = plt.subplots(figsize=(10, 10))

filterRegex = "Fil"
directory, m, a = getDirectoryMA( filterRegex=filterRegex)

posixDict = getPosixForQuantities(m, a, directory=directory)
solutionFamilyArray = openPosixDict(posixDict)
solutionFamily = {key:value[-1] for key, value in solutionFamilyArray.items()}

compute = Vacuum_Polarization(bcs=directory.name, subtractMasslessPertVacuumPolarization=False)
compute.setConfigFromDict(solutionFamily)

ooo = []
aaa = []

for n in range(1, 20):
    if n == 0:
        continue
    color = next(compute.colorCycle)
    fundamentalPert = compute.perturbativePhi(n)
    
    aaa.append(max(solutionFamily["eigenstates"][compute.maxN + n - 1])/max(fundamentalPert(compute.z)))

    # ax1.plot(solutionFamily["eigenstates"][compute.maxN + n - 1], color=color)
    # ax1.plot(fundamentalPert(compute.z), '--', color=color)


ax1.plot(aaa)
plt.show()
