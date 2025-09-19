from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt

vp = Vacuum_Polarization(
<<<<<<< HEAD
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
=======
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
>>>>>>> 56f278ba38b8e6dca00bf7b5f466caed29774e7e
        hadamard=True,
        relaxParameter=1,
        threads=5,
        )

vp.fullScript()
