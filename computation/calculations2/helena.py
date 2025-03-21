import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from Vacuum_Polarization import Vacuum_Polarization

vp = Vacuum_Polarization(lambdaMin=1.1, bisectionTol=1e-8, 
        subtractPertModes=False,
        pointsPerMode=8, 
        smoothing=False, 
        maxN=20,
        hadamard=False,
        parallelization=True,
        )

vp.stop = vp.maxN - 1

vp.singleIteration()


plt.plot(vp.z, vp.rho)
plt.plot(vp.z, vp.perturbativeModeRho(vp.stop + 1, vp.z))


plt.show()
