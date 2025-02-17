from Vacuum_Polarization import Vacuum_Polarization

import matplotlib.pyplot as plt
import sys
import numpy as np

lambdaMin = 10
lambdaMax = 26
parallelization = True
lambdaStep = 1
lambdaStepMin = 1e-8
nIterations = None
smoothing = True
showPlots = False
plotForEachLambda = False
savePlots = False
saveData = True
antisymmetric = False
a = 1

bcsArray = ["dirichlet", "neumann"]
massArray = [0, 1, 5, 10]
maxNArray = [50]


for bcs in bcsArray:
    for m in massArray:
        for maxN in maxNArray:
            directory = f"Hadamard"
            if bcs == "neumann": 
                if m == 0:
                    continue
                lambdaMin = 0.1

            computation = Vacuum_Polarization(
                    maxN=int(maxN),
                    m=m,
                    lambdaMin=lambdaMin,
                    a=1,
                    lambdaMax=lambdaMax,
                    parallelization=parallelization,
                    lambdaStep=lambdaStep,
                    lambdaStepMin=lambdaStepMin,
                    nIterations=nIterations,
                    smoothing=smoothing,
                    showPlots=showPlots,
                    plotForEachLambda=plotForEachLambda,
                    savePlots=savePlots,
                    saveData=saveData,
                    directory=directory,
                    bcs=bcs,
                    antisymmetric=antisymmetric,
                    )

            try:
                computation.fullScript()
            except Exception as e:
                raise (e)
                continue
