from Vacuum_Polarization import Vacuum_Polarization

import matplotlib.pyplot as plt
import sys
import numpy as np

m= 0

for bcs in ["dirichlet"]:
    for maxN in [50]:
        directory = f"HadamardModeSubtractionMaxN{maxN}"
        # directory = f"savingSoluitiosTest"
        computation = Vacuum_Polarization(
                maxN=int(maxN),
                m=m,
                lambdaMin = 10,
                lambdaMax =26,
                lambdaStep = 2,
                lambdaStepMin=1e-8,
                nIterations= None,
                smoothing=True,
                showPlots=False,
                plotForEachLambda=False,
                savePlots=False,
                saveData=True,
                directory=directory,
                bcs=bcs,
                )

        try:
            computation.fullScript()
        except Exception as e:
            print(e)
            continue
