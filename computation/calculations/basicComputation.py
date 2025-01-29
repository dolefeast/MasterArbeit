from Vacuum_Polarization import Vacuum_Polarization

import matplotlib.pyplot as plt
import sys
import numpy as np


m= 1

for bcs in ["dirichlet"]:
    for maxN in [33]:
        directory = f"NoFilterMaxN{maxN}lambda10"
        computation = Vacuum_Polarization(
                maxN=int(maxN),
                m=m,
                lambdaMin = 0.1,
                lambdaMax =10.1,
                lambdaStep = 10000000,
                lambdaStepMin=1e-6,
                relaxParameter=1,
                nIterations=1,
                walkbackFactor=0.2,
                iterateTol=1e-2,
                maxNIterations=300,
                a = 1, 
                smoothing=True,
                showPlots=True,
                plotForEachLambda=True,
                savePlots=False,
                saveData=False,
                directory=directory,
                extrapolateA0=True,
                bcs=bcs,
                dynamicRelaxParameter=False,
                )

        try:
            computation.fullScript()
        except Exception as e:
            raise e
