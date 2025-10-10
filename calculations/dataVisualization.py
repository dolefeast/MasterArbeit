from scipy.interpolate import CubicSpline
import sys
import os
from itertools import cycle, count
import matplotlib.pyplot as plt
import numpy as np


import matplotlib_format
from Vacuum_Polarization import Vacuum_Polarization

vp = Vacuum_Polarization()

try:
    filterRegex = sys.argv[1]
except IndexError:
    filterRegex = ""

CMAP = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
BASE_DIR = "/home/santi/Projects/Vacuum Polarization/MasterArbeit/calculations"

def checkBrokenLambdaIndex(solutionFamilyArray):
    """
    Checks for the last lambda value in which the number of modes remains constant
    """
    try:
        maxN = len(solutionFamilyArray["eigenvalues"][0])
    except TypeError:
        return None
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

def plotDesiredLambda(solutionFamilyArray, axesDict, desiredLambdaValues, lineStyle):
    colorCycle = cycle("#247ba0, #a40e4c, #28502e, #1a1b41, #8d0801".split(", "))

    lambdaValueArray = solutionFamilyArray["lambdaValue"]
    desiredLambdaIndices = [ getIndexForValueFromBelow(lambdaValueArray, value) for value in desiredLambdaValues ]

    for i in desiredLambdaIndices:
        lineColor = next(colorCycle)
        for quantity in desiredQuantities:
            axesDict[quantity].set_xlabel("z")
            axesDict[quantity].set_xlabel(r"$\lambda$")
            axesDict[quantity].set_xlabel(r"$\lambda$")
            if i is None:
                continue
            try:
                A0Induced = solutionFamilyArray[quantity][i]
            except:
                continue

            lambdaValue = lambdaValueArray[i]

            nPoints = len(A0Induced)
            z = np.linspace(0, 1, nPoints)

            axesDict[quantity].plot(
                    z,
                    A0Induced,
                    color=lineColor,
                    linestyle=lineStyle,
                    #label=f"$\lambda = {round(lambdaValue, 3)} $",
                    )

def formatAxes(axRho, ax_A0Induced, axQScreened, axQInduced, axeigenstate, axOmega, solutionFamilyArray):
    axRho, ax_A0Induced, axQScreened, axQInduced, axeigenstate, axOmega = axes
    axRho.set_xlabel(r"$z$")
    ax_A0Induced.set_xlabel(r"$z$")

    axQScreened.set_xlabel(r"$\lambda$")
    axQInduced.set_xlabel(r"$\lambda$")

    axRho.set_ylabel(r"$\rho$")
    ax_A0Induced.set_ylabel(r"$A_0^\text{br}$")
    axeigenstate.set_ylabel(r"$\phi_1$")

    axOmega.set_xlabel(r"$\lambda$")
    axOmega.set_ylabel(r"$\omega_n$")

    axQInduced.set_ylabel(r"$Q_I$")
    axQScreened.set_ylabel(r"$\lambda + Q_I$")

    try:
        maxN = len(solutionFamilyArray["eigenvalues"][0])
    except TypeError:
        maxN = 1
    if maxN == 1:
        omegaYlim = 1.05*solutionFamilyArray["eigenvalues"][0][:]
        omegaYlim[0] = 0
    else:
        omegaYlim = (0, 1.05*solutionFamilyArray["eigenvalues"][0][2])

    axOmega.set_ylim(omegaYlim)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()
    fig7.tight_layout()

def plotInducedCharge(solutionFamilyArray, axQInduced, axQScreened):
    brokenLambdaIndex = checkBrokenLambdaIndex(solutionFamilyArray)

    lambdaValueArray = solutionFamilyArray["lambdaValue"]

    inducedQList = []

    for i, lambdaValue in enumerate(solutionFamilyArray["lambdaValue"][:brokenLambdaIndex]):
        A0Induced = solutionFamilyArray["A0Induced"][i]
        inducedQList.append(maxDerivative(A0Induced))
    
    inducedTotalChargeList = [x+y for x, y in zip(lambdaValueArray[:brokenLambdaIndex], inducedQList)]

    label = directory.parent.name
    if len(lambdaValueArray)==1:
        marker = "x"
    else:
        marker = "-"

    axQInduced.plot(lambdaValueArray[:brokenLambdaIndex], inducedQList, label=label, color=color, marker=marker)
    axQScreened.plot(lambdaValueArray[:brokenLambdaIndex], inducedTotalChargeList, label=label, marker=marker, color=color)

def main(
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
    
    color = next(CMAP)
    directory, m, a, solutionFamilyArray = vp.openSolutionFamilyArray(filterRegex)

    plotInducedCharge(solutionFamilyArray, axQInduced, axQScreened)

    lambdaValueArray = solutionFamilyArray["lambdaValue"]
        
    print("The maximum lambda value obtained is", max(lambdaValueArray[:brokenLambdaIndex])) #, ", with a minimum lambda step of", lambdaValueArray[-1]-lambdaValueArray[-2])

    axOmega.plot(lambdaValueArray[:brokenLambdaIndex], solutionFamilyArray["eigenvalues"][:brokenLambdaIndex], color=color,)
    axOmega.plot([], [],  label=label, color=color,)

    plotDesiredLambda(solutionFamilyArray, axesDict, desiredLambdaValues, lineStyle)
    
    formatAxes(axRho, ax_A0Induced, axQScreened, axQInduced, axeigenstate, axOmega, solutionFamilyArray)

    return solutionFamilyArray

maxLambdaDensity = None  # There are more solutions for one region of lambda than there are for other, this approximately keeps the density of solutions per lambda value fixed for visualization purposes

fig1 = plt.figure(f'vacuum polarization')  #,  dpi=600)
fig2 = plt.figure(f'induced potential')  #,  dpi=600)
fig3 = plt.figure(f'mode energy')  #,  dpi=600)
fig4 = plt.figure(f'Screening charge')  #,  dpi=600)
fig5 = plt.figure(f'Screened charge')  #,  dpi=600)
fig6 = plt.figure(f'Eigenstates')  #,  dpi=600)
fig7 = plt.figure(f'Mode charge density')  #,  dpi=600)

axRho = fig1.subplots()
ax_A0Induced = fig2.subplots()
axOmega = fig3.subplots()
axQInduced = fig4.subplots()
axQScreened = fig5.subplots()
axeigenstate = fig6.subplots()
axModeChargeDensity = fig7.subplots()

desiredQuantities = ["rho", "A0Induced", "eigenstates"]

axesDict = {"rho":axRho, "A0Induced":ax_A0Induced, "eigenstates":axeigenstate}

desiredLambdaValues = [1,]

formats = cycle(["-", "--", ":", "-."])

for nCurves in count():
    try:
        solutionFamilyArray = main(
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

if nCurves > 1:
    axQScreened.legend(loc="best")
    axQInduced.legend(loc="best")
    axOmega.legend(loc="best")
    ncol = 2
else:
    ncol = 1

axRho.legend(loc="best", ncol=ncol, markerscale=0.3)
ax_A0Induced.legend(loc="best", ncol=ncol, markerscale=0.3)
axeigenstate.legend(loc="best", ncol=ncol, markerscale=0.3)

date = time.strftime("%Y%m%d-%H%M%S")

if not os.path.exists(date):
    os.makedirs(f"{BASE_DIR}/figures/{date}")

fig1.savefig(f"{BASE_DIR}/figures/{date}/rhoEvolution.pdf")
fig2.savefig(f"{BASE_DIR}/figures/{date}/A0InducedEvolution.pdf")
fig3.savefig(f"{BASE_DIR}/figures/{date}/eigenvalues.pdf")
fig4.savefig(f"{BASE_DIR}/figures/{date}/inducedCharge.pdf")
fig5.savefig(f"{BASE_DIR}/figures/{date}/electricFieldInduced.pdf")
fig6.savefig(f"{BASE_DIR}/figures/{date}/eigenstate.pdf")
fig7.savefig(f"{BASE_DIR}/figures/{date}/modeChargeDensity.pdf")