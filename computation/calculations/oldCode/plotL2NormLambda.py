from plotScripts.readDataScripts import getPosixForQuantities, openPosixDict, strToFloat
from plotScripts.readFiles import getDirectoryMA

import matplotlib.pyplot as plt

bcs = "neumann"
maxLambdaDensity=None

directory, m, a = getDirectoryMA(bcs)

desiredQuantities = ["rho", "A0Induced"]
posixDict = getPosixForQuantities(m, a, directory=directory, maxLambdaDensity=maxLambdaDensity)

solutionFamilyArray = openPosixDict(posixDict, desiredQuantities=desiredQuantities)

lambdaValueArray = solutionFamilyArray["lambdaValue"]

L2Norms = {
        "rho": [], 
        "A0Induced": [], 
        }

fig = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)}')
(axRho, axA0Induced) = fig.subplots(2)
axesDict = {"rho":axRho, "A0Induced":axA0Induced}

for quantity in desiredQuantities:
    axesDict[quantity].set_xlabel("$\lambda$")
    for i, lambdaValue in enumerate(lambdaValueArray):
            L2Norms[quantity].append(
                    # rootMeanSquare(
                    # solutionFamilyArray[quantity][i]
                    # )
                    max(solutionFamilyArray[quantity][i])
                    )
    axesDict[quantity].plot(lambdaValueArray, L2Norms[quantity])

axRho.set_ylabel(r"max($\rho$)")
axA0Induced.set_ylabel(r"max($-\int \int \rho$)")
# axRho.setYlabel(r"$\sqrt{\langle \rho^2\rangle}$")
# ax_A0Induced.setYlabel(r"$\sqrt{\langle A_0^2\rangle}$")

fig.tight_layout()
plt.show()
