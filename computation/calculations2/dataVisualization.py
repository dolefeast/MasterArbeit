from scipy.interpolate import CubicSpline
import sys
import os
from itertools import cycle, count
import matplotlib.pyplot as plt
import numpy as np
import time

import scienceplots

plt.style.use(["science", "high-contrast"])
plt.rcParams["figure.figsize"] = (3.5, 2.6)
plt.rcParams["font.size"] = "5.4"
plt.rcParams["axes.labelsize"] = "13"
plt.rcParams["xtick.labelsize"] = "13"
plt.rcParams["ytick.labelsize"] = "13"
plt.rcParams["lines.linewidth"] = "0.9"

from Vacuum_Polarization import Vacuum_Polarization

vp = Vacuum_Polarization()

try:
    filterRegex = sys.argv[1]
except IndexError:
    filterRegex = ""

cmap = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])

def checkBrokenLambdaIndex(solutionFamilyArray):
    """
    Checks for the last lambda value in which the number of modes remains constant
    """
    maxN = len(solutionFamilyArray["eigenvalues"][0])
    for brokenLambdaIndex, lambdaValue in enumerate(solutionFamilyArray["eigenvalues"]):
        tempMaxN = len(solutionFamilyArray["eigenvalues"][brokenLambdaIndex])
        if maxN != tempMaxN:
            return brokenLambdaIndex-1

    return None

def getIndexForValueFromBelow(array, value):
    """
    Given an array, return the biggest index i for which array[i] <= value
    """
    for i, element in enumerate(array):
        if element > value:
            return i-1
    return -1

def maxDerivative(y):
    x = np.linspace(0, 1, len(y))
    y = CubicSpline(x, y)

    return -y.derivative()(1/2)

def plotOneSolutionFamily(
        lineStyle,
        desiredLambdaValues,
        filterRegex, 
        axRho,
        ax_A0Induced,
        axOmega,
        axQInduced,
        axQScreened,
        axesDict,
            ):

    colorCycle = cycle("#247ba0, #a40e4c, #28502e, #1a1b41, #8d0801".split(", "))
    directory, m, a, solutionFamilyArray = vp.openSolutionFamilyArray(filterRegex)

    inducedQList = []

    maxN = len(solutionFamilyArray["eigenvalues"][0])
    if maxN == 1: 
        omegaYlim = 1.05*solutionFamilyArray["eigenvalues"][0][:]
        omegaYlim[0] = 0
    else:
        omegaYlim = (0, 1.05*solutionFamilyArray["eigenvalues"][0][2])

    omegaYlim = (0, 10)

    axOmega.set_ylim(omegaYlim)
    lambdaValueArray = solutionFamilyArray["lambdaValue"]
        
    brokenLambdaIndex = checkBrokenLambdaIndex(solutionFamilyArray)

    color = next(cmap)

    for i, lambdaValue in enumerate(solutionFamilyArray["lambdaValue"][:brokenLambdaIndex]):
        recoveredQuantity = solutionFamilyArray["A0Induced"][i]
        inducedQList.append(maxDerivative(recoveredQuantity))
    
    inducedTotalChargeList = [x+ y for x, y in zip(lambdaValueArray[:brokenLambdaIndex], inducedQList)]

    label = f"{str(directory.parent.name)}"
    axQInduced.plot(lambdaValueArray[:brokenLambdaIndex], inducedQList, label=label, color=color)
    axQScreened.plot(lambdaValueArray[:brokenLambdaIndex], inducedTotalChargeList, label=label, color=color)

    print("The maximum lambda value obtained is", max(lambdaValueArray[:brokenLambdaIndex])) #, ", with a minimum lambda step of", lambdaValueArray[-1]-lambdaValueArray[-2])

    axOmega.plot(lambdaValueArray[:brokenLambdaIndex], solutionFamilyArray["eigenvalues"][:brokenLambdaIndex], color=color, label=label)

    desiredLambdaIndex = [ getIndexForValueFromBelow(lambdaValueArray, value) for value in desiredLambdaValues ]
    print(desiredLambdaIndex)

    for i in desiredLambdaIndex:
        lineColor = next(colorCycle)
        for quantity in desiredQuantities:
            axesDict[quantity].set_xlabel("z")
            axesDict[quantity].set_xlabel(r"$\lambda$")
            axesDict[quantity].set_xlabel(r"$\lambda$")
            if i is None:
                continue
            recoveredQuantity = solutionFamilyArray[quantity][i]

            lambdaValue = lambdaValueArray[i]

            nPoints = len(recoveredQuantity)
            z = np.linspace(0, 1, nPoints) 

            axesDict[quantity].plot(
                    z,
                    recoveredQuantity,
                    # color=lineColor,
                    # label=f"$\lambda={lambdaValueArray[i]}$", 
                    label=f"$\lambda = {lambdaValue} $", 
                    )


    axRho.set_xlabel(r"$z$")
    ax_A0Induced.set_xlabel(r"$z$")

    axQScreened.set_xlabel(r"$\lambda$")
    axQInduced.set_xlabel(r"$\lambda$")

    axRho.set_ylabel(r"$\rho$")
    ax_A0Induced.set_ylabel(r"$A_0^\text{br}$")

    axOmega.set_xlabel(r"$\lambda$")
    axOmega.set_ylabel(r"$\omega_n$")

    axQInduced.set_ylabel(r"$Q_I$")
    axQScreened.set_ylabel(r"$\lambda + Q_I$")

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()

    return solutionFamilyArray

maxLambdaDensity = None  # There are more solutions for one region of lambda than there are for other, this approximately keeps the density of solutions per lambda value fixed for visualization purposes


fig1 = plt.figure(f'vacuum polarization')  #,  dpi=600)
fig2 = plt.figure(f'induced potential')  #,  dpi=600)
fig3 = plt.figure(f'mode energy')  #,  dpi=600)
fig4 = plt.figure(f'Screening charge')  #,  dpi=600)
fig5 = plt.figure(f'Screened charge')  #,  dpi=600)


axRho = fig1.subplots()
ax_A0Induced = fig2.subplots()
axOmega = fig3.subplots()
axQInduced = fig4.subplots()
axQScreened = fig5.subplots()

desiredQuantities = ["rho", "A0Induced"]

axesDict = {"rho":axRho, "A0Induced":ax_A0Induced}

desiredLambdaValues = [14.4, 16.01, 18.01]


formats = cycle(["-"])#, "--", ":", "-."])

for n in count():
    try:
        solutionFamilyArray = plotOneSolutionFamily(
            next(formats),
            desiredLambdaValues,
            filterRegex, 
            axRho,
            ax_A0Induced,
            axOmega,
            axQInduced,
            axQScreened,
            axesDict,
            )
       

    except TypeError as e:
        print(e)
        break


if n > 1:
    axQScreened.legend(loc="best")
    axQInduced.legend(loc="best")
    axOmega.legend(loc="best")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    axOmega.legend(by_label.values(), by_label.keys())
    ncol = 2
else:
    ncol = 1

axRho.legend(loc="best", ncol=ncol, markerscale=0.3)
ax_A0Induced.legend(loc="best", ncol=ncol, markerscale=0.3)

date = time.strftime("%Y%m%d-%H%M%S")


plt.show()

if not os.path.exists(date):
    os.makedirs("figures/" + date)

fig1.savefig(f"figures/{date}/vacuumPolarizationEvolution.pdf")
fig2.savefig(f"figures/{date}/A0InducedEvolution.pdf")
fig3.savefig(f"figures/{date}/eigenvalues.pdf")
fig4.savefig(f"figures/{date}/inducedCharge.pdf")
fig5.savefig(f"figures/{date}/electricFieldInduced.pdf")

