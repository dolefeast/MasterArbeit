from Vacuum_Polarization import Vacuum_Polarization

import matplotlib.pyplot as plt
import sys
import numpy as np

lambdaMax = 30
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
maxNIterations = 50
a = 1

bcsArray = ["dirichlet"]
massArray = [5, 10]
maxNArray = [25, 50, 75, 100, 125]

for bcs in bcsArray:
    for m in massArray:
        for maxN in maxNArray:
            directory = f"Hadamard"
            if bcs == "neumann": 
                if m == 0:
                    continue

            lambdaMin = 1
            computation = Vacuum_Polarization(
                    maxN=int(maxN),
                    m=m,
                    lambdaMin=lambdaMin,
                    a=a,
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
                    maxNIterations=maxNIterations,
                    )

            try:
                computation.fullScript()
            except Exception as e:
                print(e)
                continue
