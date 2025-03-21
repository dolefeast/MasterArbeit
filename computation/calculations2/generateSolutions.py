import matplotlib.pyplot as plt
import numpy as np

from Vacuum_Polarization import Vacuum_Polarization

fig1, axRho = plt.subplots()
fig2, axA0 = plt.subplots()
fig3, axmaxA0 = plt.subplots()

pendiente = []
modeArray = []

modeArray = list(range(180, 10, -20))

vp = Vacuum_Polarization(lambdaMin=.1, bisectionTol=1e-8, 
        subtractPertModes=False,
        pointsPerMode=8, 
        smoothing=True, 
        hadamard=True,
        maxN=modeArray[0],
        )

vp.calculateEigenstatesParallel()

for maxN in modeArray:

    vp.maxN = maxN

    vp.eigenstates = vp.eigenstates[:maxN]
    vp.eigenvalues = vp.eigenvalues[:maxN]

    vp.calculateRho()
    if vp.smoothing:
        vp.convolveRho()
    vp.calculateA0Induced()

    vp.directory = f"Filtering! N {vp.maxN}"
    vp.saveDataScript()

