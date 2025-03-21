from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt

vp = Vacuum_Polarization(parallelization=True, bisectionTol=1e-8, lambdaMin=10, a=1, nIterations=None, maxN=60, saveData=True, directory="Self-consistent solutions", showPlots=False, smoothing=True)

vp.fullScript()
