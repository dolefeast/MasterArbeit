from Vacuum_Polarization import Vacuum_Polarization
from plotScripts.readFiles import getDirectoryMA
from plotScripts.readDataScripts import getPosixForQuantities, openPosixDict, strToFloat

import matplotlib.pyplot as plt


fig1, ax1 = plt.subplots(figsize=(10, 10))

def antisymmetrize(arr:list):
    return [(xForward - xBackwards)/2 for xForward, xBackwards in zip(arr, arr[::-1])]

for _ in range(1):
    filterRegex = "NoFi*lam"
    directory, m, a = getDirectoryMA( filterRegex=filterRegex)
    posixDict = getPosixForQuantities(m, a, directory=directory)
    solutionFamilyArray = openPosixDict(posixDict)
    solutionFamily = {key:value[-1] for key, value in solutionFamilyArray.items()}
    compute = Vacuum_Polarization(bcs=directory.name)

    compute.setConfigFromDict(solutionFamily)
    x, filteredRho = compute.extendAndFilter(compute.rho)

    # plt.plot(compute.rho, label=r"unfiltered $\rho$")
    ax1.plot(compute.z, filteredRho,  label=r"filtered $\rho$, maxN="+str(compute.maxN))
    # plt.plot(antisymmetrize(filteredRho), label="antisymmetrized signal")

plt.legend(loc="best")

fig1.suptitle(f"$\lambda = {compute.lambdaValue}$")
ax1.set_xlabel("z")
plt.show()
