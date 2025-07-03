from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt

vp = Vacuum_Polarization(
        parallelization=True,
        bisectionTol=1e-8,
        lambdaMin=0.5,
        lambdaStep=0.5,
        lambdaMax=49999,
        a=0,
        nIterations=None,
        maxN=5,
        saveData=True,
        directory="External field approximation Neumann",
        showPlots=False,
        smoothing=True,
        m=1,
        bcs="neumann",
        hadamard=True,
        relaxParameter=1,
        threads=5,
        )

vp.fullScript()
