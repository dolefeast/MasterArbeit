import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from Vacuum_Polarization import Vacuum_Polarization

vp = Vacuum_Polarization(maxN = 1, lambdaMin=1.123189379, bisectionTol=1e-5, 
        subtractPertModes=True,
        pointsPerMode=8, 
        smoothing=False, 
        hadamard=True,
        saveData=True,
        directory="Testing",
        )

vp.singleIteration()

vp.saveDataScript()
