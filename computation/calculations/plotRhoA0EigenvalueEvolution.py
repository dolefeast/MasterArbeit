from plotScripts.readDataScripts import openSolutionFamilyArray

import sys
from itertools import cycle
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

cmap = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

def plotOneSolutionFamily(
        filterRegex, 
        axRho,
        ax_A0Induced,
        axOmega,
        axL2NormsRho,
        axL2NormsA0Induced,
            ):
    directory, solutionFamilyArray = openSolutionFamilyArray(filterRegex)

    L2NormsDict = {"rho":[], "A0Induced":[]}

    maxN = len(solutionFamilyArray["eigenvalues"][0])//2
    if maxN == 1: 
        omegaYlim = 1.2*solutionFamilyArray["eigenvalues"][0][:]
        omegaYlim[0] = 0
    else:
        omegaYlim = (0, 10)

    axOmega.set_ylim(omegaYlim)
    lambdaValueArray = solutionFamilyArray["lambdaValue"]
        
    maxN = len(solutionFamilyArray["eigenvalues"][0])

    for brokenLambdaIndex, lambdaValue in enumerate(lambdaValueArray):
        tempMaxN = len(solutionFamilyArray["eigenvalues"][brokenLambdaIndex])
        if maxN != tempMaxN:
            break
    maxN = tempMaxN

    color = next(cmap)

    for quantity in desiredQuantities:
        for i, lambdaValue in enumerate(solutionFamilyArray["lambdaValue"][:brokenLambdaIndex]):
            recoveredQuantity = solutionFamilyArray[quantity][i]
            L2NormsDict[quantity].append(max(recoveredQuantity))
        axesL2Dict[quantity].plot(lambdaValueArray[:brokenLambdaIndex], L2NormsDict[quantity], label=f"{str(directory.parent.name)}/{str(directory.name)}", color=color)
        axesL2Dict[quantity].legend(loc="best")
        axesL2Dict[quantity].set_xlabel(r"$\lambda$")

    print("The maximum lambda value obtained is", lambdaValueArray[brokenLambdaIndex]) #, ", with a minimum lambda step of", lambdaValueArray[-1]-lambdaValueArray[-2])

    axOmega.plot(lambdaValueArray[:brokenLambdaIndex], solutionFamilyArray["eigenvalues"][:brokenLambdaIndex], color=color, label = f"{str(directory.parent.name)}/{str(directory.name)}")

    return solutionFamilyArray



maxLambdaDensity = None

# brokenLambdaIndex -= 2


fig1 = plt.figure(f'vacuum polarization', figsize=(6, 4))
fig2 = plt.figure(f'induced potential', figsize=(6, 4))
fig3 = plt.figure(f'mode energy', figsize=(6, 4))
fig4 = plt.figure(f'L2 polarization', figsize=(6, 4))
fig5 = plt.figure(f'L2 potential', figsize=(6, 4))

axRho = fig1.subplots()
ax_A0Induced = fig2.subplots()
axOmega = fig3.subplots()
axL2NormsRho = fig4.subplots()
axL2NormsA0Induced = fig5.subplots()

desiredQuantities = ["rho", "A0Induced"]

axesDict = {"rho":axRho, "A0Induced":ax_A0Induced}
axesL2Dict = {"rho":axL2NormsRho, "A0Induced":axL2NormsA0Induced}

while True:
    try:
        solutionFamilyArray = plotOneSolutionFamily(
            filterRegex, 
            axRho,
            ax_A0Induced,
            axOmega,
            axL2NormsRho,
            axL2NormsA0Induced)
       
        for quantity in desiredQuantities:
            axesDict[quantity].set_xlabel("z")
            axesDict[quantity].set_xlabel(r"$\lambda$")
            for i, lambdaValue in enumerate(solutionFamilyArray["lambdaValue"][-2:]):
                recoveredQuantity = solutionFamilyArray[quantity][i]

                nPoints = len(recoveredQuantity)
                z = np.linspace(0, 1, nPoints) 

                axesDict[quantity].plot(
                        z,
                        recoveredQuantity,
                        label=f"$\lambda={lambdaValue}$", 
                        )

    except TypeError as e:
        break




axRho.set_xlabel(r"$z$", fontsize=17)
ax_A0Induced.set_xlabel(r"$z$", fontsize=17)

axRho.set_ylabel(r"$\rho$", fontsize=17)
ax_A0Induced.set_ylabel(r"$A_0$", fontsize=17)

axOmega.set_xlabel(r"$\lambda$", fontsize=17)
axOmega.set_ylabel(r"$\omega_N$", fontsize=17)
axOmega.legend(loc="best")


# fig1.tight_layout()
# fig2.tight_layout()
# fig3.tight_layout()
# fig4.tight_layout()

axL2NormsRho.set_ylabel(r"$L_2(\rho)$")

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axOmega.legend(by_label.values(), by_label.keys())

plt.show()
