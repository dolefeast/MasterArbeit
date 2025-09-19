from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt

vp = Vacuum_Polarization(
        bisectionTol=1e-8,
        lambdaMin=0.5,
        lambdaStep=1,
        lambdaMax=25,
        walkbackFactor=0.5,
        a=1,
        nIterations=None,
        maxN=15,
        saveData=True,
        directory="External field approximation",
        showPlots=False,
        smoothing=True,
        m=0,
        bcs="dirichlet",
        hadamard=True,
        relaxParameter=1,
        threads=5,
        )

vp.fullScript()
