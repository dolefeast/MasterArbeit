from plotScripts.readDataScripts import getPosixForQuantities, openPosixDict, strToFloat
from plotScripts.readFiles import getDirectoryMA

import sys
import matplotlib.pyplot as plt
import numpy as np

def getIndexForValueFromBelow(array, value):
    """
    Given an array, return the biggest index i for which array[i] <= value
    """
    for i, element in enumerate(array):
        if element > value:
            return i-1

try:
    filterRegex = sys.argv[1]
except IndexError:
    filterRegex = ""


bcs = "neumann"
maxLambdaDensity = None

directory, m, a = getDirectoryMA(bcs, filterRegex=filterRegex)

desiredQuantities = ["eigenstates"]
posixDict = getPosixForQuantities(m, a, directory=directory, maxLambdaDensity=maxLambdaDensity, bcs=bcs)
solutionFamilyArray = openPosixDict(posixDict, desiredQuantities=["eigenvalues"] + desiredQuantities)

maxN = len(solutionFamilyArray["eigenstates"][0])//2

lambdaValueArray = solutionFamilyArray["lambdaValue"]

print("The maximum lambda value obtained is", lambdaValueArray[-1]) #, ", with a minimum lambda step of", lambdaValueArray[-1]-lambdaValueArray[-2])

fig1 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)}1')

axRho = fig1.subplots()

axesDict = {"eigenstates":axRho}

nSolutions = len(lambdaValueArray )

for quantity in desiredQuantities:
    axesDict[quantity].set_xlabel("z")
    for i, lambdaValue in enumerate(lambdaValueArray):
        recoveredQuantity = solutionFamilyArray[quantity][i][maxN]

        nPoints = len(recoveredQuantity)
        z = np.linspace(0, 1, nPoints) 

        axesDict[quantity].plot(
                z,
                recoveredQuantity,
                label=f"$\lambda={lambdaValue}$", 
                )


axRho.set_ylabel(r"$\rho$", fontsize=17)

fig1.tight_layout()

plt.show()
