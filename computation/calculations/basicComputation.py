from Vacuum_Polarization import Vacuum_Polarization

import matplotlib.pyplot as plt
import sys
import numpy as np


m= 0

for bcs in ["dirichlet"]:
    for maxN in [75]:
        directory = f"HadamardModeSubtractionMaxN{maxN}"
        computation = Vacuum_Polarization(
                maxN=int(maxN),
                m=m,
                lambdaMin = 10,
                lambdaMax =26,
                lambdaStep = 2,
                lambdaStepMin=1e-8,
                relaxParameter=1,
                nIterations=None,
                walkbackFactor=0.2,
                iterateTol=1e-2,
                maxNIterations=300,
                a = 1, 
                smoothing=True,
                showPlots=False,
                plotForEachLambda=False,
                savePlots=False,
                saveData=True,
                directory=directory,
                extrapolateA0=True,
                bcs=bcs,
                dynamicRelaxParameter=False,
                subtractMasslessPertVacuumPolarization=True
                )

        try:
            computation.fullScript()
        except Exception as e:
            raise e
