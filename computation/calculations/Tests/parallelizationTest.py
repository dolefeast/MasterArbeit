# Testing calculating the eigenstates through parallelization

from time import time
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

from Vacuum_Polarization import Vacuum_Polarization
from differentRhoComparison import calculateNorm

def calculateEigens(compute:object, omegaLowerUpper:tuple)->float:
    if compute.bcs == "neumann":
        initialValues = (1, 0)
        bcsIndex = 1
    elif compute.bcs == "dirichlet":
        initialValues = (0, 1)
        bcsIndex = 0

    parameterizedODE = lambda omega: solve_ivp(lambda z, y: compute.KleinGordonEquation(z, y, omega),
            t_span=(0, 1),
            y0=initialValues, 
            dense_output=True)

    try:
        omega = compute.findRootBisection(lambda omega: parameterizedODE(omega).sol(1)[bcsIndex], *omegaLowerUpper)
    except RuntimeError:
        return None, None
    
    normCalc = calculateNorm(omega, parameterizedODE(omega), compute.lambdaValue)

    return  parameterizedODE(omega).sol(compute.z)[0]/np.sqrt(abs(normCalc))

def ccccccccccccc(omegaLowerUpper):
    return calculateEigens(compute, omegaLowerUpper)

compute = Vacuum_Polarization(maxN = 20)

omegaUpperArray = compute.bisectionMethodUpperBound(compute.eigenvalues)
omegaLowerArray = compute.bisectionMethodLowerBound(compute.eigenvalues)

t0 = time()
with Pool() as pool:
    eigenstates = pool.map(ccccccccccccc, zip(omegaLowerArray, omegaUpperArray))

t1 = time()

compute.calculateEigenstates()
compute.normalizeEigenstates()
t2 = time()

print(f"""Comparing the eigenstate family runtime for parallel vs serial for maxN={compute.maxN}
            Using multiprocessing.pool: t={t1-t0}s
            Using serial calculation: t={t2-t1}s"""
            )

plt.plot(eigenstates[compute.maxN], label="parallelization")
plt.plot(compute.eigenstates[compute.maxN], label="serial")

plt.show()
