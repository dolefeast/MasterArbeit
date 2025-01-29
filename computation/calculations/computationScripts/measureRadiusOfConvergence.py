from Vacuum_Polarization import Vacuum_Polarization
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt  

import os
import sys

def getAxLims(ax):
    return (ax.getXlim(), ax.getYlim())

def fixedPointPlot(fun):

    def inner():

        A0InducedValue = compute.A0Induced(1)
        A0FullValue = compute.lambdaValue/2 + compute.A0(1)
        fun()
        fixedPointLineInduced.append((A0InducedValue, compute.A0Induced(1)))
        fixedPointLineFull.append((A0FullValue, compute.lambdaValue/2 + compute.A0(1)))
        # ax.plot(a1, compute.A0Induced(1), 'o', color=compute.color, label=f"$A0Induced(1)={compute.A0Induced(1)}$")

    return inner

compute = Vacuum_Polarization(
        maxN=int(sys.argv[1]),
        iterateTol = 0.0005,
        lambdaMin=10,
        lambdaMax=17,
        lambdaStep=6,
        save=False,
        savePlots=False,
        showPlots=False,
        smoothing=True,
        ambjorn=False,
        m=0,
        relaxParameter=0.6,
        maxNiterations=500,
        )    

# I know lambdaValue = 16 is a stable solution
# want to measure the radius of convergence
# So: Get to lambdaValue = 16, store the configuration
# of the system in that case, and then study what happens
# when using different starting A0_induced

try:
    compute.fullScript()
except KeyboardInterrupt:
    pass

# It should have stopped at a configuration that yields a bifurcation
# because it breaks in case of RecursionError

bifurcationConfig = dict()
bifurcationConfig["eigenvalues"] = compute.eigenvalues
bifurcationConfig["eigenstates"] = compute.eigenstates
bifurcationConfig["A0Induced"] = compute.A0Induced(compute.z)

bifurcationConfig["lambdaValue"] = 17

# Looking for the radius of convergence of a certain configuration
# This means: What's the furthest away I can start from the converged
# A0 so that my problem converges?

relaxParameterArray = [0.001, 0.01, 0.05, 0.1][::-1]
A0FactorArray = np.arange(-15, 15, 2)

for c in relaxParameterArray:
    print(f"c = {c}")

    successArray = []
    fig1, (axA0InducedFixedPoint, axA0InducedIteration) = plt.subplots(nrows=1, ncols=2, figsize=(16,9))
    fig2, (axA0FullFixedPoint, axA0FullIteration) = plt.subplots(nrows=1, ncols=2, figsize=(16,9))
    fig1.suptitle(f"c = {c}")
    fig2.suptitle(f"c = {c}")

    compute.relaxParameter = c
    compute.iterateTol = c/10
    compute.coreIteration = fixedPointPlot(compute.coreIteration)

    for A0Factor in A0FactorArray:
        directory = f"convergence_radiusMaxn_{compute.maxn}"

        fixedPointLineInduced = []
        fixedPointLineFull = []

        compute.setConfigFromDict(bifurcationConfig)
        compute.A0InducedHistory = [
                [
                CubicSpline(
                compute.z, 
                A0Factor * bifurcationConfig["A0Induced"]
                )
                    ]
                ] # Need to overwrite this

        try:
            compute.color = next(compute.colorCycle)
            compute.constantLambdaIterations()
            success = compute.A0Induced(1)
        except KeyboardInterrupt:
            break
        except RuntimeError as e: 

            print(e)
            success = False
        finally:
            if len(fixedPointLineInduced) != 0:
                axA0InducedIteration.plot([x for x, y in fixedPointLineInduced], 
                        'x-',
                        label=f"A0 induced(1) = {round(A0Factor * bifurcationConfig['A0Induced'][-1], 3)}",
                        color=compute.color)

                axA0InducedFixedPoint.plot(
                        [x  for x, y in fixedPointLineInduced],
                        [y  for x, y in fixedPointLineInduced], 
                        'x-',
                        label=f"A0 induced(1) = {round(A0Factor * bifurcationConfig['A0Induced'][-1], 3)}",
                        color=compute.color)

                axA0FullIteration.plot([x for x, y in fixedPointLineFull], 
                        'x-',
                        label=f"A0 full(1) = {round(fixedPointLineFull[0][0], 3)}",
                        color=compute.color)

                axA0FullFixedPoint.plot(
                        [x  for x, y in fixedPointLineFull],
                        [y  for x, y in fixedPointLineFull], 
                        'x-',
                        label=f"A0 full(1) = {round(fixedPointLineFull[0][0], 3)}",
                        color=compute.color)

        
        successArray.append(success)
    axA0InducedIteration.set_xlabel("k")
    axA0InducedIteration.set_ylabel("$A_k(1)$ induced")

    axA0FullIteration.set_xlabel("k")
    axA0FullIteration.set_ylabel(f"full $A_k(1)$ + {compute.lambdaValue/2}")

    axA0InducedFixedPoint.set_xlabel("$A_k(1)$ induced")
    axA0InducedFixedPoint.set_ylabel("$A_{k+1}(1)$ induced")

    axA0FullFixedPoint.set_xlabel("full $A_{{k}}(1)$ + {}".format(compute.lambdaValue/2))
    axA0FullFixedPoint.set_ylabel("full $A_{{k+1}}(1)$ + {}".format(compute.lambdaValue/2))

    axA0InducedIteration.legend(loc="best")
    axA0FullIteration.legend(loc="best")
    axA0InducedFixedPoint.legend(loc="best")
    axA0FullFixedPoint.legend(loc="best")

    axA0InducedIteration.legend(loc="best")
    axA0FullIteration.legend(loc="best")
    axA0InducedFixedPoint.legend(loc="best")
    axA0FullFixedPoint.legend(loc="best")

    axFixedPointList = [axA0InducedFixedPoint, axA0FullFixedPoint]

    for ax in axFixedPointList:
        axLims = getAxLims(ax) 
        print(axLims)
        ax.plot([-999, -999], [999, 999], 'g--', alpha=0.2)
        ax.set_xlim(axLims[0])
        ax.set_ylim(axLims[1])

    try:
        os.mkdir(f"figures/{directory}")
    except FileExistsError:
        pass

    fig1.savefig(f"figures/{directory}/fixed_point_A0_induced_c_{compute.float_to_str(c)}.pdf")
    fig2.savefig(f"figures/{directory}/fixed_point_A0_full_c_{compute.float_to_str(c)}.pdf")
