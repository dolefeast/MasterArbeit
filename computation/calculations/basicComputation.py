from Vacuum_Polarization import Vacuum_Polarization
from computationScripts.eigenstatesPlot import main

import matplotlib.pyplot as plt
import sys
import numpy as np

directory = "noFilterMaxN12"

main(
        maxN=12,
        m=0,
        lambdaMin = 1,
        lambdaMax = 12,
        lambdaStep = 2,
        lambdaStepMin=1e-8,
        relaxParameter=1,
        walkbackFactor=0.2,
        iterateTol=1e-2,
        nIterations=1,
        maxNIterations=300,
        a = 0, 
        smoothing=True,
        showPlots=True,
        plotForEachLambda=True,
        savePlots=True,
        saveData=False,
        directory=directory,
        extrapolateA0=False,
        bcs="dirichlet",
        dynamicRelaxParameter=False,
        )
