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


maxLambdaDensity = None

directory, m, a = getDirectoryMA(filterRegex=filterRegex)

desiredQuantities = ["rho", "A0Induced"]
posixDict = getPosixForQuantities(m, a, directory=directory, maxLambdaDensity=maxLambdaDensity)
solutionFamilyArray = openPosixDict(posixDict, desiredQuantities=["eigenvalues"] + desiredQuantities)

maxN = len(solutionFamilyArray["eigenvalues"][0])//2
if maxN == 1: 
    omegaYlim = 1.2*solutionFamilyArray["eigenvalues"][0][:]
    omegaYlim[0] = 0
else:
    omegaYlim = (0, 10)

lambdaValueArray = solutionFamilyArray["lambdaValue"]

print("The maximum lambda value obtained is", lambdaValueArray[-1]) #, ", with a minimum lambda step of", lambdaValueArray[-1]-lambdaValueArray[-2])

fig1 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)} vacuum polarization')
fig2 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)} induced potential')
fig3 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)} mode energy')
fig4 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)} L2 polarization')
fig5 = plt.figure(f'{directory}, m={strToFloat(m)}, a={strToFloat(a)} L2 potential')

axRho = fig1.subplots()
ax_A0Induced = fig2.subplots()
axOmega = fig3.subplots()
axL2NormsRho = fig4.subplots()
axL2NormsA0Induced = fig5.subplots()

L2NormsDict = {"rho":[], "A0Induced":[]}
axesDict = {"rho":axRho, "A0Induced":ax_A0Induced}
axesL2Dict = {"rho":axL2NormsRho, "A0Induced":axL2NormsA0Induced}

for quantity in desiredQuantities:
    axesDict[quantity].set_xlabel("z")
    axesDict[quantity].set_xlabel(r"$\lambda$")
    for i, lambdaValue in enumerate(lambdaValueArray):
        recoveredQuantity = solutionFamilyArray[quantity][i]

        nPoints = len(recoveredQuantity)
        z = np.linspace(0, 1, nPoints) 

        axesDict[quantity].plot(
                z,
                recoveredQuantity,
                label=f"$\lambda={lambdaValue}$", 
                )
        L2NormsDict[quantity].append(max(recoveredQuantity))
    axesL2Dict[quantity].plot(lambdaValueArray, L2NormsDict[quantity])

axOmega.plot(lambdaValueArray, solutionFamilyArray["eigenvalues"], 'b')

axRho.set_ylabel(r"$\rho$", fontsize=17)
ax_A0Induced.set_ylabel(r"$A_0$", fontsize=17)

axOmega.set_xlabel(r"$\lambda$", fontsize=17)
axOmega.set_ylabel(r"$\omega_N$", fontsize=17)

axOmega.set_ylim(omegaYlim)

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()

plt.show()
